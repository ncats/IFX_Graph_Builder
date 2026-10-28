import os
import pickle
import re
import csv
import gzip
import json
from decimal import Decimal
from types import SimpleNamespace
from typing import List, Dict, Generator

from sqlalchemy import or_

from src.constants import Prefix, CHEMBL_PATENT_SOURCE_ID, CHEMBL_FUNCTIONAL_ASSAY_CODE, CHEMBL_BINDING_ASSAY_CODE, \
    CHEMBL_SMALL_MOLECULE_CODE, CHEMBL_SINGLE_PROTEIN_CODE, HUMAN_TAX_ID, DataSourceName
from src.input_adapters.chembl.tables import Activities, CompoundRecords, CompoundStructures, MoleculeDictionary, \
    Assays, TargetDictionary, ComponentSequence, TargetComponents, Docs
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.registry_dataset import RegistryDatasetKind
from src.input_adapters.sql_adapter import MySqlAdapter
from src.interfaces.input_adapter import InputAdapter
from src.models.ligand import Ligand, ProteinLigandEdge, ActivityDetails
from src.models.node import EquivalentId
from src.models.protein import Protein
from src.shared.db_credentials import DBCredentials


class ChemblAdapter(MySqlAdapter):
    pchembl_cutoff: float
    version_info: DatasourceVersionInfo

    def __init__(self, data_source, pchembl_cutoff: float = 5, credentials: DBCredentials = None):
        self.data_source = data_source
        self.version_info = data_source.version_info()
        self.pchembl_cutoff = pchembl_cutoff
        if data_source.kind is RegistryDatasetKind.EXTERNAL:
            if credentials is None:
                raise ValueError("External ChEMBL database input requires credentials")
            MySqlAdapter.__init__(self, credentials)
        elif data_source.kind is RegistryDatasetKind.SOURCE:
            if pchembl_cutoff < 5:
                raise ValueError(
                    "ChEMBL activity export contains only pChEMBL >= 5; "
                    f"configured cutoff is {pchembl_cutoff}"
                )
            self.snapshot_cutoff = Decimal(str(pchembl_cutoff))
            contract = json.loads(
                data_source.file("chembl_activity_export_contract.json").read_text(
                    encoding="utf-8"
                )
            )
            if (contract.get("snapshot_id") != data_source.snapshot_id
                    or contract.get("pchembl_cutoff") != 5
                    or contract.get("null_encoding")
                    != "backslash_N_with_leading_backslash_escape"
                    or not isinstance(contract.get("row_count"), int)
                    or contract["row_count"] <= 0):
                raise ValueError("ChEMBL activity export contract does not match its snapshot")
            self.expected_row_count = contract["row_count"]
        else:
            raise ValueError(f"Unsupported ChEMBL dataset kind: {data_source.kind}")

    def _activity_cache_file(self, pchembl_cutoff: float) -> str:
        schema = self.credentials.schema or "default"
        version = self.version_info.version or "unknown"
        safe_schema = re.sub(r"[^A-Za-z0-9_.-]+", "_", schema)
        safe_version = re.sub(r"[^A-Za-z0-9_.-]+", "_", version)
        return f"input_files/chembl/activities_{safe_schema}_{safe_version}_{pchembl_cutoff}.pkl"

    def fetch_activity_data(self, pchembl_cutoff: float):
        if self.data_source.kind is RegistryDatasetKind.SOURCE:
            return self._snapshot_activity_rows()
        cache_file = self._activity_cache_file(pchembl_cutoff)
        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as f:
                print('loading chembl data from cache')
                return pickle.load(f)

        query = self.get_all_activities_query(pchembl_cutoff)

        results = []
        for row in query.yield_per(50000):
            results.append(row)

        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, 'wb') as f:
            pickle.dump(results, f)

        return results

    def _snapshot_activity_rows(self):
        path = self.data_source.file("chembl_activities_pchembl5.tsv.gz")
        required = {
            "activity_id", "compound_id", "compound_name", "canonical_smiles",
            "target_id", "uniprot_id", "action_type", "pchembl_value",
            "standard_type", "assay_type", "pubmed_id", "src_id",
        }
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(
                    f"ChEMBL activity export {path} lacks required columns: "
                    + ", ".join(sorted(missing))
                )
            row_count = 0
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(
                        f"Malformed ChEMBL activity export row in {path} at line {line_number}"
                    )
                try:
                    row = {
                        key: (None if value == "\\N" else value[1:]
                              if value.startswith("\\\\") else value)
                        for key, value in row.items()
                    }
                    if not row["compound_id"] or not row["target_id"]:
                        raise ValueError("blank compound or target ID")
                    value = Decimal(row["pchembl_value"])
                    if value < 5:
                        raise ValueError("pChEMBL value below export cutoff")
                    activity_id = int(row["activity_id"])
                    pubmed_id = int(row["pubmed_id"]) if row["pubmed_id"] else None
                    src_id = int(row["src_id"])
                    row_count += 1
                    if value < self.snapshot_cutoff:
                        continue
                    yield SimpleNamespace(
                        activity_id=activity_id,
                        compound_id=row["compound_id"],
                        compound_name=row["compound_name"],
                        canonical_smiles=row["canonical_smiles"],
                        target_id=row["target_id"],
                        uniprot_id=row["uniprot_id"],
                        action_type=row["action_type"],
                        pchembl_value=value,
                        standard_type=row["standard_type"],
                        assay_type=row["assay_type"],
                        pubmed_id=pubmed_id,
                        src_id=src_id,
                    )
                except (ValueError, TypeError, ArithmeticError) as error:
                    raise ValueError(
                        f"Invalid ChEMBL activity export row in {path} "
                        f"at line {line_number}: {error}"
                    ) from error
            if row_count != self.expected_row_count:
                raise ValueError(
                    f"ChEMBL activity export {path} has {row_count} rows; "
                    f"contract declares {self.expected_row_count}"
                )

    def get_all_activities_query(self, pchembl_cutoff):
        return (self.get_session().query(
            Activities.activity_id,
            MoleculeDictionary.chembl_id.label('compound_id'),
            CompoundRecords.compound_name,
            CompoundStructures.canonical_smiles,
            TargetDictionary.chembl_id.label('target_id'),
            ComponentSequence.accession.label('uniprot_id'),
            TargetDictionary.pref_name,
            Activities.action_type,
            Activities.pchembl_value,
            Activities.standard_type,
            Assays.assay_type,
            Docs.journal,
            Docs.year,
            Docs.volume,
            Docs.issue,
            Docs.first_page,
            Docs.pubmed_id,
            CompoundRecords.src_id
        ).join(CompoundRecords, Activities.molregno == CompoundRecords.molregno)
                .join(CompoundStructures, Activities.molregno == CompoundStructures.molregno)
                .join(MoleculeDictionary, Activities.molregno == MoleculeDictionary.molregno)
                .join(Assays, Activities.assay_id == Assays.assay_id)
                .join(TargetDictionary, Assays.tid == TargetDictionary.tid)
                .join(TargetComponents, TargetDictionary.tid == TargetComponents.tid)
                .join(ComponentSequence, TargetComponents.component_id == ComponentSequence.component_id)
                .join(Docs, Activities.doc_id == Docs.doc_id, isouter=True)
                .filter(
            MoleculeDictionary.structure_type == CHEMBL_SMALL_MOLECULE_CODE,
            TargetDictionary.target_type == CHEMBL_SINGLE_PROTEIN_CODE,
            ComponentSequence.tax_id == HUMAN_TAX_ID,
        )
                .filter(
            Activities.standard_flag == 1,
            Activities.pchembl_value >= pchembl_cutoff,
            Activities.standard_relation == '=',
            Assays.assay_type.in_([CHEMBL_FUNCTIONAL_ASSAY_CODE, CHEMBL_BINDING_ASSAY_CODE]),
            or_(
                Activities.doc_id.isnot(None),
                CompoundRecords.src_id == CHEMBL_PATENT_SOURCE_ID
            )
        )
                .distinct()
                )


class DrugNodeAdapter(InputAdapter, ChemblAdapter):

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.ChEMBL

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[List[Ligand], None, None]:
        activity_results = self.fetch_activity_data(self.pchembl_cutoff)

        drug_dict: Dict[str, Ligand] = {}

        for row in activity_results:
            id = EquivalentId(id=row.compound_id, type=Prefix.CHEMBL_COMPOUND)
            id_str = id.id_str()

            if id_str not in drug_dict:
                drug_dict[id_str] = Ligand(
                    id=id.id_str(),
                    smiles=row.canonical_smiles,
                    name=row.compound_name
                )

        yield list(drug_dict.values())


class ProteinDrugEdgeAdapter(InputAdapter, ChemblAdapter):
    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.ChEMBL

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[List[ProteinLigandEdge], None, None]:
        activity_results = self.fetch_activity_data(self.pchembl_cutoff)

        relationship_map: Dict[str, ProteinLigandEdge] = {}

        for row in activity_results:
            compound_id = EquivalentId(id=row.compound_id, type=Prefix.CHEMBL_COMPOUND)
            uniprot = row.uniprot_id
            if uniprot is not None and len(uniprot) > 0:
                pro_id = EquivalentId(id=uniprot, type=Prefix.UniProtKB)
            else:
                target_id = row.target_id
                pro_id = EquivalentId(id=target_id, type=Prefix.CHEMBL_PROTEIN)

            rel_id = f"{compound_id.id_str()}|{pro_id.id_str()}"
            if rel_id in relationship_map:
                pro_lig_edge = relationship_map[rel_id]
            else:
                lig = Ligand(
                    id=compound_id.id_str()
                )
                pro = Protein(id = pro_id.id_str())
                pro_lig_edge = ProteinLigandEdge(
                    start_node=pro, end_node=lig
                )
                relationship_map[rel_id] = pro_lig_edge

            activity_details = ActivityDetails(
                ref_id=row.activity_id,
                activity_source=DataSourceName.ChEMBL,
                act_value=float(row.pchembl_value),
                act_type=row.standard_type,
                action_type=row.action_type,
                assay_type=row.assay_type,
                act_pmids=[row.pubmed_id],
                comment=f"Chembl Source ID: {row.src_id}"
            )
            pro_lig_edge.details.append(activity_details)

        yield list(relationship_map.values())
