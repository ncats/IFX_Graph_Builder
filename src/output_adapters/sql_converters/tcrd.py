import hashlib
from datetime import date, datetime
from typing import Union, List, Optional
from src.constants import Prefix, DataSourceName
from src.models.drgc_resource import DRGCResource
from src.models.disease import Disease, DiseaseParentEdge, MondoDiseaseParentEdge, DODiseaseParentEdge, ProteinDiseaseEdge, TINXImportanceEdge
from src.models.dto_class import DTOClass, DTOClassParentEdge, ProteinDTOClassEdge
from src.models.expression import ProteinTissueExpressionEdge
from src.models.external_link import ExternalLinkProvider, ProteinExternalLinkEdge
from src.models.generif import GeneGeneRifEdge
from src.models.go_term import GoType, GoTerm, GoTermHasParent, ProteinGoTermEdge
from src.models.harmonizome import HarmonizomeGeneAttributeType, HarmonizomeHgramCDF
from src.models.keyword import ProteinKeywordEdge
from src.models.ligand import Ligand, ProteinLigandEdge
from src.models.node import EquivalentId
from src.models.ncats_datasource import NcatsDataSourceInfo, NcatsDataSourceMapEntry
from src.models.panther_class import PantherClass, ProteinPantherClassEdge
from src.models.pathway import ProteinPathwayEdge
from src.models.ppi import PPIEdge
from src.models.protein import Protein
from src.models.gwas_trait import ProteinGwasTraitEdge, GwasTrait
from src.models.mouse_phenotype import OrthologGeneMousePhenotypeEdge, ProteinMousePhenotypeEdge
from src.models.tcrd_disease_ontology import MondoTerm, MondoTermParentEdge, DOTerm, DOTermParentEdge
from src.models.tissue import Tissue, TissueParentEdge
from src.models.ortholog import OrthologGene, ProteinOrthologGeneEdge
from src.models.virus import Virus, ViralPPIDetail, ViralPPIEdge, ViralProtein
from src.models.word_count import WordCount
from src.input_adapters.shared.hcop import HCOP_TAXID_TO_SPECIES
from src.shared.sqlalchemy_tables.pharos_tables_new import (
    DrgcResource as mysqlDrgcResource, Protein as mysqlProtein, Xref, Alias, Target, TDL_info, T2TC, GO, GOParent, GoA,
    GeneRif, GeneRif2Pubmed, Protein2Pubmed, Ligand as mysqlLigand, LigandActivity,
    Uberon, UberonParent, Tissue as mysqlTissue, Expression, Gtex,
    Mondo, MondoParent, MondoXref, Disease as mysqlDisease, DiseaseType, DO, DOParent,
    NcatsDisease, NcatsD2DA, Pathway as mysqlPathway, PantherClass as mysqlPantherClass, P2PC, PPI as mysqlPPI,
    Tiga as mysqlTiga, TigaProvenance, TinxImportance,
    DTO as mysqlDTO, DTOParent, P2DTO, Pmscore, NhProtein as mysqlNhProtein, Phenotype as mysqlPhenotype,
    Ortholog as mysqlOrtholog, PatentCount, Ptscore, Virus as mysqlVirus, ViralProtein as mysqlViralProtein,
    ViralPPI as mysqlViralPPI, Affiliate, ExtLink,
    GeneAttributeType as mysqlGeneAttributeType, HgramCDF as mysqlHgramCDF,
    NcatsDataSource, NcatsDataSourceMap,
    WordCount as mysqlWordCount,
)
from src.output_adapters.sql_converters.output_converter_base import SQLOutputConverter
from src.shared.sqlalchemy_tables.pharos_tables_new import Base as TCRDBase


class TCRDOutputConverter(SQLOutputConverter):

    def __init__(self):
        super().__init__(sql_base=TCRDBase)
        self._known_disease_types: set = set()
        self._disease_name_by_id: dict[str, str] = {}
        self._known_mondo_ids: set[str] = set()
        self._disease_mondoid_by_id: dict[str, str | None] = {}
        self._next_ncats_disease_id: int | None = None
        self._seen_protein2pubmed: set[tuple[int, str, int | None, str]] = set()
        self._target_id_by_uniprot: dict[str, int] = {}
        self._warned_missing_drgc_targets: set[tuple[str, str | None]] = set()
        self._protein_id_by_uniprot: dict[str, int] = {}
        self._protein_id_by_geneid: dict[int, int] = {}
        self._ambiguous_protein_geneids: set[int] = set()
        self._protein_id_by_symbol: dict[str, int] = {}
        self._ambiguous_protein_symbols: set[str] = set()
        self._warned_missing_harmonizome_proteins: set[tuple[int | None, str | None, int | None, str | None]] = set()
        self._converters = {
            # Protein
            Protein: [self.protein_converter, self.target_converter, self.t2tc_converter,
                      self.protein_alias_converter, self.protein_xref_converter, self.tdl_info_converter,
                      self.pmscore_converter, self.ptscore_converter, self.patent_count_converter, self.generif_from_publications_converter,
                      self.generif_assoc_from_publications_converter, self.p2p_from_publications_converter],
            # GeneRif
            GeneGeneRifEdge: [self.generif_converter, self.generif_assoc_converter, self.p2p_converter],
            # GO
            GoTerm: [self.goterm_converter],
            GoTermHasParent: [self.goterm_parent_converter],
            ProteinGoTermEdge: [self.goa_converter],
            # Ligand
            Ligand: [self.ligand_converter],
            ProteinLigandEdge: [self.ligand_edge_converter],
            # Tissue / Expression
            Tissue: [self.uberon_converter],
            TissueParentEdge: [self.uberon_parent_converter],
            ProteinTissueExpressionEdge: [self.tissue_lookup_converter,
                                          self.expression_converter,
                                          self.gtex_converter],
            # Direct MySQL side-loads
            DRGCResource: [self.drgc_resource_converter],
            HarmonizomeGeneAttributeType: [self.harmonizome_gene_attribute_type_converter],
            HarmonizomeHgramCDF: [self.harmonizome_hgram_cdf_converter],
            OrthologGene: [self.nhprotein_converter],
            ProteinOrthologGeneEdge: [self.ortholog_converter],
            OrthologGeneMousePhenotypeEdge: [self.phenotype_converter],
            ProteinMousePhenotypeEdge: [self.protein_mouse_phenotype_converter],
            MondoTerm: [self.mondo_table_converter, self.mondo_xref_converter],
            MondoTermParentEdge: [self.mondo_parent_table_converter],
            DOTerm: [self.do_table_converter],
            DOTermParentEdge: [self.do_parent_table_converter],
            # Disease
            Disease: [self.ncats_disease_converter],
            DiseaseParentEdge: None,
            MondoDiseaseParentEdge: None,
            DODiseaseParentEdge: None,
            ProteinDiseaseEdge: [self.disease_type_converter, self.disease_converter, self.ncats_d2da_converter],
            TINXImportanceEdge: [self.tinx_importance_converter],
            # TIGA
            GwasTrait: None,
            ProteinGwasTraitEdge: [self.tiga_converter, self.tiga_provenance_converter],
            # Pathway
            ProteinPathwayEdge: [self.pathway_converter],
            # PPI
            PPIEdge: [self.ppi_converter],
            # Viral PPI
            Virus: [self.virus_converter],
            ViralProtein: [self.viral_protein_converter],
            ViralPPIEdge: [self.viral_ppi_converter],
            # Panther
            PantherClass: [self.panther_class_converter],
            ProteinPantherClassEdge: [self.p2pc_converter],
            # DTO
            DTOClass: [self.dto_converter],
            DTOClassParentEdge: [self.dto_parent_converter],
            ProteinDTOClassEdge: [self.p2dto_converter],
            # Keyword
            ProteinKeywordEdge: [self.keyword_xref_converter],
            # Simple linkouts
            ExternalLinkProvider: [self.affiliate_converter],
            ProteinExternalLinkEdge: [self.extlink_converter],
            # Publication word-cloud background
            WordCount: [self.word_count_converter],
            # Pharos compatibility datasource tables
            NcatsDataSourceInfo: [self.ncats_datasource_converter],
            NcatsDataSourceMapEntry: [self.ncats_datasource_map_converter],
        }

    def get_object_converters(self, obj_cls) -> Union[callable, List[callable], None]:
        return self._converters.get(obj_cls)

    def resolve_id(self, table, id):
        if table != 'ncats_disease':
            return super().resolve_id(table, id)
        mapping = self.id_mapping.setdefault(table, {})
        if id not in mapping:
            if self._next_ncats_disease_id is None:
                self._next_ncats_disease_id = max(mapping.values(), default=0) + 1
            mapping[id] = self._next_ncats_disease_id
            self._next_ncats_disease_id += 1
        return mapping[id]

    def word_count_converter(self, obj: dict) -> mysqlWordCount:
        return mysqlWordCount(word=obj["word"], count=obj["count"])

    def ncats_datasource_converter(self, obj: dict) -> NcatsDataSource:
        return NcatsDataSource(
            dataSource=obj["dataSource"],
            dataSourceDescription=obj.get("dataSourceDescription"),
            url=obj.get("url"),
            license=obj.get("license"),
            licenseURL=obj.get("licenseURL"),
            citation=obj.get("citation"),
        )

    def ncats_datasource_map_converter(self, obj: dict) -> NcatsDataSourceMap:
        protein_id = obj.get("protein_id")
        if protein_id is None and obj.get("protein_ifx_id"):
            protein_id = self.resolve_id("protein", obj["protein_ifx_id"])

        ligand_id = obj.get("ncats_ligand_id")
        if ligand_id is None and obj.get("ligand_identifier"):
            ligand_id = self.resolve_id("ligand", obj["ligand_identifier"])

        return NcatsDataSourceMap(
            dataSource=obj["dataSource"],
            url=obj.get("url"),
            license=obj.get("license"),
            licenseURL=obj.get("licenseURL"),
            protein_id=protein_id,
            ncats_ligand_id=ligand_id,
            disease_name=obj.get("disease_name"),
        )

    @staticmethod
    def _parse_patent_family_token(token: str) -> tuple[int, int] | None:
        if not isinstance(token, str) or ":" not in token:
            return None
        year_str, family_id_str = token.split(":", 1)
        try:
            return int(year_str), int(family_id_str)
        except ValueError:
            return None

    def get_preload_queries(self, session):
        return [{
            "table": 'protein',
            "data": session.query(mysqlProtein.id, mysqlProtein.ifx_id).all()
        }, {
            "table": 'ncats_disease_mondoid',
            "data": session.query(NcatsDisease.id, NcatsDisease.mondoid)
            .filter(NcatsDisease.mondoid.isnot(None))
            .all()
        }, {
            "table": 'tissue',
            "data": session.query(mysqlTissue.id, mysqlTissue.name).all()
        }, {
            "table": 'panther_class',
            "data": session.query(mysqlPantherClass.id, mysqlPantherClass.pcid).all()
        }, {
            "table": 'protein_by_uniprot',
            "data": session.query(mysqlProtein.id, mysqlProtein.uniprot).filter(mysqlProtein.uniprot.isnot(None)).all()
        }, {
            "table": 'target_by_uniprot',
            "data": session.query(Target.id, mysqlProtein.uniprot)
            .join(T2TC, T2TC.target_id == Target.id)
            .join(mysqlProtein, mysqlProtein.id == T2TC.protein_id)
            .all()
        }, {
            "table": 'ligand',
            "data": session.query(mysqlLigand.id, mysqlLigand.identifier).all()
        }]

    def preload_id_mappings(self, session):
        super().preload_id_mappings(session)
        self._target_id_by_uniprot = dict(self.id_mapping.get("target_by_uniprot", {}))
        self._protein_id_by_uniprot = dict(self.id_mapping.get("protein_by_uniprot", {}))
        self._protein_id_by_geneid, self._ambiguous_protein_geneids = self._build_unambiguous_geneid_map(session)
        self._protein_id_by_symbol, self._ambiguous_protein_symbols = self._build_unambiguous_symbol_map(session)
        self._known_disease_types = {
            row[0] for row in session.query(DiseaseType.name).all() if row[0]
        }
        self._known_mondo_ids = {
            row[0] for row in session.query(Mondo.mondoid).all() if row[0]
        }

    @staticmethod
    def _build_unambiguous_geneid_map(session) -> tuple[dict[int, int], set[int]]:
        geneid_rows = (
            session.query(mysqlProtein.geneid, mysqlProtein.id)
            .filter(mysqlProtein.geneid.isnot(None))
            .all()
        )
        geneid_to_id: dict[int, int] = {}
        ambiguous: set[int] = set()
        for geneid, protein_id in geneid_rows:
            if geneid is None:
                continue
            key = int(geneid)
            if key in ambiguous:
                continue
            if key in geneid_to_id and geneid_to_id[key] != protein_id:
                ambiguous.add(key)
                del geneid_to_id[key]
                continue
            geneid_to_id[key] = protein_id
        return geneid_to_id, ambiguous

    @staticmethod
    def _build_unambiguous_symbol_map(session) -> tuple[dict[str, int], set[str]]:
        symbol_rows = (
            session.query(mysqlProtein.sym, mysqlProtein.id)
            .filter(mysqlProtein.sym.isnot(None))
            .all()
        )
        symbol_to_id: dict[str, int] = {}
        ambiguous: set[str] = set()
        for symbol, protein_id in symbol_rows:
            if not symbol:
                continue
            key = symbol.upper()
            if key in ambiguous:
                continue
            if key in symbol_to_id and symbol_to_id[key] != protein_id:
                ambiguous.add(key)
                del symbol_to_id[key]
                continue
            symbol_to_id[key] = protein_id
        return symbol_to_id, ambiguous

    def _resolve_legacy_harmonizome_protein_id(self, obj: dict) -> int | None:
        uniprot = obj.get("legacy_uniprot")
        if uniprot and uniprot in self._protein_id_by_uniprot:
            return self._protein_id_by_uniprot[uniprot]

        geneid = obj.get("legacy_geneid")
        if geneid is not None:
            try:
                geneid_key = int(geneid)
            except (TypeError, ValueError):
                geneid_key = None
            if (
                geneid_key is not None
                and geneid_key not in self._ambiguous_protein_geneids
                and geneid_key in self._protein_id_by_geneid
            ):
                return self._protein_id_by_geneid[geneid_key]

        symbol = obj.get("legacy_symbol")
        if symbol:
            symbol_key = symbol.upper()
            if symbol_key not in self._ambiguous_protein_symbols and symbol_key in self._protein_id_by_symbol:
                return self._protein_id_by_symbol[symbol_key]

        warning_key = (obj.get("legacy_protein_id"), uniprot, geneid, symbol)
        if warning_key not in self._warned_missing_harmonizome_proteins:
            self._warned_missing_harmonizome_proteins.add(warning_key)
            print(
                "Skipping Harmonizome row with no destination protein match "
                f"(legacy protein_id={obj.get('legacy_protein_id')}, "
                f"uniprot={uniprot}, geneid={geneid}, symbol={symbol})"
            )
        return None

    # --- Protein ---

    def protein_converter(self, obj: dict) -> mysqlProtein:
        gene_ids = [p for p in obj['xref'] if p.startswith('NCBIGene:')]
        gene_id = gene_ids[0].split(':')[-1] if gene_ids else None
        ensembl_ids = [p for p in obj['xref'] if p.startswith('ENSEMBL:ENSP')]
        string_id = ensembl_ids[0].split(':')[-1] if ensembl_ids and gene_id is None else None
        protein_id = self.resolve_id('protein', obj['id'])
        uniprot_id = obj['uniprot_id']
        self._target_id_by_uniprot[uniprot_id] = protein_id
        self._protein_id_by_uniprot[uniprot_id] = protein_id
        if gene_id is not None:
            geneid_key = int(gene_id)
            existing_id = self._protein_id_by_geneid.get(geneid_key)
            if existing_id is not None and existing_id != protein_id:
                self._ambiguous_protein_geneids.add(geneid_key)
                del self._protein_id_by_geneid[geneid_key]
            elif geneid_key not in self._ambiguous_protein_geneids:
                self._protein_id_by_geneid[geneid_key] = protein_id
        if obj.get('symbol'):
            symbol_key = obj['symbol'].upper()
            existing_id = self._protein_id_by_symbol.get(symbol_key)
            if existing_id is not None and existing_id != protein_id:
                self._ambiguous_protein_symbols.add(symbol_key)
                del self._protein_id_by_symbol[symbol_key]
            elif symbol_key not in self._ambiguous_protein_symbols:
                self._protein_id_by_symbol[symbol_key] = protein_id
        return mysqlProtein(
            id=protein_id,
            ifx_id=obj['id'],
            name=obj.get('gene_name'),
            description=obj.get('name'),
            uniprot=uniprot_id,
            sym=obj.get('symbol'),
            geneid=gene_id,
            stringid=string_id,
            seq=obj['sequence'],
            dtoid=(obj.get('dtoid') or '').replace(':', '_') or None,
            dtoclass=obj.get('dtoclass'),
            provenance=obj['provenance'],
            preferred_symbol=obj.get('preferred_symbol'),
            novelty=min(obj.get('novelty') or []) if obj.get('novelty') else None,
        )

    # --- Harmonizome compatibility tables ---

    @staticmethod
    def harmonizome_gene_attribute_type_converter(obj: dict) -> mysqlGeneAttributeType:
        return mysqlGeneAttributeType(
            id=obj["legacy_id"],
            name=obj["name"],
            association=obj["association"],
            description=obj["description"],
            resource_group=obj["resource_group"],
            measurement=obj["measurement"],
            attribute_group=obj["attribute_group"],
            attribute_type=obj["attribute_type"],
            pubmed_ids=obj.get("pubmed_ids"),
            url=obj.get("url"),
            provenance=obj["provenance"],
        )

    def harmonizome_hgram_cdf_converter(self, obj: dict) -> mysqlHgramCDF | None:
        protein_id = self._resolve_legacy_harmonizome_protein_id(obj)
        if protein_id is None:
            return None

        return mysqlHgramCDF(
            protein_id=protein_id,
            type=obj["type"],
            attr_count=obj["attr_count"],
            attr_cdf=obj["attr_cdf"],
            provenance=obj["provenance"],
        )

    def tdl_info_converter(self, obj: dict) -> List[TDL_info]:
        tdl_infos = []
        uniprot_description = obj.get('uniprot_function')
        if uniprot_description:
            tdl_infos.append(TDL_info(
                itype="UniProt Function",
                protein_id=self.resolve_id('protein', obj['id']),
                string_value=uniprot_description
            ))
        ncbi_gene_summary = obj.get('ncbi_gene_summary')
        if ncbi_gene_summary:
            tdl_infos.append(TDL_info(
                itype="NCBI Gene Summary",
                protein_id=self.resolve_id('protein', obj['id']),
                string_value=ncbi_gene_summary
            ))
        if obj.get('antibody_count') and len(obj['antibody_count']) > 0:
            antibody_count = max([int(p) for p in obj['antibody_count']])
            if antibody_count > 0:
                tdl_infos.append(TDL_info(
                    itype="Ab Count",
                    protein_id=self.resolve_id('protein', obj['id']),
                    integer_value=antibody_count
                ))
        if obj.get('pm_score') and len(obj['pm_score']) > 0:
            pm_score = max([float(p) for p in obj['pm_score']])
            tdl_infos.append(TDL_info(
                itype="JensenLab PubMed Score",
                protein_id=self.resolve_id('protein', obj['id']),
                number_value=pm_score
            ))
        if obj.get('pt_score') and len(obj['pt_score']) > 0:
            pt_score = max([float(p) for p in obj['pt_score']])
            tdl_infos.append(TDL_info(
                itype="PubTator Score",
                protein_id=self.resolve_id('protein', obj['id']),
                number_value=pt_score
            ))
        if obj.get('patent_family_mentions') and len(obj['patent_family_mentions']) > 0:
            unique_families = {
                parsed[1]
                for parsed in (self._parse_patent_family_token(token) for token in obj['patent_family_mentions'])
                if parsed is not None
            }
            if unique_families:
                tdl_infos.append(TDL_info(
                    itype="SureChEMBL Patent Family Count",
                    protein_id=self.resolve_id('protein', obj['id']),
                    integer_value=len(unique_families)
                ))
        return tdl_infos

    def protein_alias_converter(self, obj: dict) -> List[Alias]:
        protein_id = self.resolve_id('protein', obj['id'])
        ids = [EquivalentId.parse(x) for x in obj['xref']]
        aliases = []
        for s in [x for x in ids if x.type in (Prefix.OldSymbol, Prefix.Symbol)]:
            aliases.append(Alias(protein_id=protein_id, type='symbol', value=s.id))
        for u in [x for x in ids if x.type == Prefix.UniProtKB]:
            aliases.append(Alias(protein_id=protein_id, type='uniprot', value=u.id))
        for n in [x for x in ids if x.type == Prefix.NCBIGene]:
            aliases.append(Alias(protein_id=protein_id, type='NCBI Gene ID', value=n.id))
        return aliases

    def protein_xref_converter(self, obj: dict) -> List[Xref]:
        protein_id = self.resolve_id('protein', obj['id'])
        ids = [EquivalentId.parse(x) for x in obj['xref']]
        return [
            Xref(
                protein_id=protein_id,
                xtype="Ensembl" if a.type == Prefix.ENSEMBL else a.type.value,
                value=a.id
            )
            for a in ids
        ]

    def pmscore_converter(self, obj: dict) -> List[Pmscore]:
        protein_id = self.resolve_id('protein', obj['id'])
        rows = []
        for entry in obj.get('pm_score_by_year') or []:
            year = entry.get('year')
            score = entry.get('score')
            if year is None or score is None:
                continue
            rows.append(
                Pmscore(
                    protein_id=protein_id,
                    year=int(year),
                    score=float(score),
                )
                )
        return rows

    def ptscore_converter(self, obj: dict) -> List[Ptscore]:
        protein_id = self.resolve_id('protein', obj['id'])
        rows = []
        for entry in obj.get('pt_score_by_year') or []:
            year = entry.get('year')
            score = entry.get('score')
            if year is None or score is None:
                continue
            rows.append(
                Ptscore(
                    protein_id=protein_id,
                    year=int(year),
                    score=float(score),
                )
            )
        return rows

    def patent_count_converter(self, obj: dict) -> List[PatentCount]:
        protein_id = self.resolve_id('protein', obj['id'])
        counts_by_year = {}
        for token in obj.get('patent_family_mentions') or []:
            parsed = self._parse_patent_family_token(token)
            if parsed is None:
                continue
            year, family_id = parsed
            counts_by_year.setdefault(year, set()).add(family_id)

        return [
            PatentCount(
                protein_id=protein_id,
                year=year,
                count=len(family_ids),
            )
            for year, family_ids in sorted(counts_by_year.items())
        ]

    def t2tc_converter(self, obj: dict) -> T2TC:
        return T2TC(
            target_id=self.resolve_id('protein', obj['id']),
            protein_id=self.resolve_id('protein', obj['id']),
            provenance=obj['provenance'])

    def target_converter(self, obj: dict) -> Target:
        return Target(
            id=self.resolve_id("protein", obj['id']),
            ifx_id=obj['id'],
            name=obj['name'],
            ttype='Single Protein',
            fam=obj['idg_family'],
            tdl=obj['tdl'],
            provenance=obj['provenance']
        )

    def drgc_resource_converter(self, obj: dict) -> mysqlDrgcResource | None:
        target_id = self._target_id_by_uniprot.get(obj["uniprot_id"])
        if target_id is None:
            warning_key = (obj["uniprot_id"], obj.get("rssid"))
            if warning_key not in self._warned_missing_drgc_targets:
                self._warned_missing_drgc_targets.add(warning_key)
                print(
                    f"Skipping DRGC resource with no destination target_id for UniProt {obj['uniprot_id']} "
                    f"(legacy target_id={obj.get('legacy_target_id')}, rssid={obj.get('rssid')})"
                )
            return None

        return mysqlDrgcResource(
            rssid=obj["rssid"],
            resource_type=obj["resource_type"],
            target_id=target_id,
            json=obj["json"],
            provenance=obj["provenance"],
        )

    # --- Ortholog / phenotype ---

    @staticmethod
    def _parse_optional_int(value) -> Optional[int]:
        if value in (None, ""):
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _extract_xref_value(obj: dict, prefix: str) -> Optional[str]:
        for xref in obj.get("xref") or []:
            if isinstance(xref, str) and xref.startswith(f"{prefix}:"):
                return xref.split(":", 1)[1]
        return None

    @staticmethod
    def _passes_legacy_ortholog_support_rule(support_sources: List[str] | None) -> bool:
        if not support_sources:
            return False
        accepted = {"Inparanoid", "OMA", "EggNOG"}
        return len(accepted.intersection(set(support_sources))) >= 2

    @staticmethod
    def _legacy_ortholog_sources_string(support_sources: List[str] | None) -> Optional[str]:
        ordered = ["Inparanoid", "OMA", "EggNOG"]
        support_set = set(support_sources or [])
        selected = [source for source in ordered if source in support_set]
        if len(selected) < 2:
            return None
        return ", ".join(selected)

    @staticmethod
    def _legacy_ortholog_mod_url(taxid: Optional[int], db_id: Optional[str], symbol: Optional[str]) -> Optional[str]:
        if taxid == 10090:
            return f"http://www.informatics.jax.org/marker/{db_id or '-'}"
        if taxid == 10116:
            rgd_id = (db_id or "").split(":", 1)[1] if db_id and ":" in db_id else "-"
            return f"http://rgd.mcw.edu/rgdweb/report/gene/main.html?id={rgd_id or '-'}"
        if taxid == 7955:
            return f"http://zfin.org/{db_id or '-'}"
        if taxid == 7227 and db_id:
            return f"http://flybase.org/reports/{db_id}.html"
        if taxid == 6239 and symbol:
            return f"http://www.wormbase.org/search/gene/{symbol}"
        if taxid == 4932 and db_id:
            return f"https://www.yeastgenome.org/locus/{db_id}"
        return None

    def nhprotein_converter(self, obj: dict) -> mysqlNhProtein:
        taxid = self._parse_optional_int(obj.get("species"))
        symbol = obj.get("symbol")
        description = obj.get("name")
        fallback_name = symbol or description or obj["id"]
        return mysqlNhProtein(
            id=self.resolve_id("nhprotein", obj["id"]),
            uniprot=self._extract_xref_value(obj, Prefix.UniProtKB.value),
            name=fallback_name,
            description=description,
            sym=symbol,
            species=HCOP_TAXID_TO_SPECIES.get(str(obj.get("species")), str(obj.get("species") or "")),
            taxid=taxid if taxid is not None else -1,
            geneid=self._parse_optional_int(obj.get("entrez_gene_id")),
            provenance=obj.get("provenance"),
        )

    def ortholog_converter(self, obj: dict) -> mysqlOrtholog | None:
        sources = self._legacy_ortholog_sources_string(obj.get("support_sources"))
        if sources is None:
            return None

        taxid = self._parse_optional_int(obj.get("species"))
        db_id = None
        source_db_ids = obj.get("source_db_ids") or []
        if source_db_ids:
            db_id = source_db_ids[0]

        symbol = None
        ortholog_symbols = obj.get("ortholog_symbols") or []
        if ortholog_symbols:
            symbol = ortholog_symbols[0]
        if not symbol:
            symbol = (obj.get("end_node") or {}).get("symbol")

        name = None
        ortholog_names = obj.get("ortholog_names") or []
        if ortholog_names:
            name = ortholog_names[0]
        if not name:
            name = (obj.get("end_node") or {}).get("name")

        species_name = HCOP_TAXID_TO_SPECIES.get(str(obj.get("species")), str(obj.get("species") or ""))

        return mysqlOrtholog(
            protein_id=self.resolve_id("protein", obj["start_id"]),
            taxid=taxid if taxid is not None else -1,
            species=species_name,
            db_id=db_id,
            geneid=self._parse_optional_int((obj.get("end_node") or {}).get("entrez_gene_id")),
            symbol=symbol or "",
            name=name or "",
            mod_url=self._legacy_ortholog_mod_url(taxid, db_id, symbol),
            sources=sources,
            provenance=obj.get("provenance"),
        )

    def phenotype_converter(self, obj: dict) -> List[mysqlPhenotype]:
        protein_links = obj.get("protein_links") or []
        term_name = (obj.get("end_node") or {}).get("name")
        rows: List[mysqlPhenotype] = []

        for detail in obj.get("details") or []:
            source = detail.get("source") or DataSourceName.IMPC.value
            if source == DataSourceName.MGI.value:
                ptype = "JAX/MGI Human Ortholog Phenotype"
                eligible_protein_ids = [
                    link.get("protein_node_id")
                    for link in protein_links
                    if link.get("protein_node_id")
                ]
                nhprotein_id = None
            else:
                ptype = DataSourceName.IMPC.value
                eligible_protein_ids = [
                    link.get("protein_node_id")
                    for link in protein_links
                    if self._passes_legacy_ortholog_support_rule(link.get("support_sources"))
                ]
                nhprotein_id = self.resolve_id("nhprotein", obj["start_id"])

            if not eligible_protein_ids:
                continue

            for protein_ifx_id in set(eligible_protein_ids):
                rows.append(
                    mysqlPhenotype(
                        ptype=ptype,
                        protein_id=self.resolve_id("protein", protein_ifx_id),
                        nhprotein_id=nhprotein_id,
                        trait=detail.get("trait"),
                        top_level_term_id=detail.get("top_level_term_id"),
                        top_level_term_name=detail.get("top_level_term_name"),
                        term_id=obj["end_id"],
                        term_name=term_name,
                        term_description=None,
                        p_value=detail.get("p_value"),
                        percentage_change=detail.get("percentage_change"),
                        effect_size=detail.get("effect_size"),
                        procedure_name=detail.get("procedure_name"),
                        parameter_name=detail.get("parameter_name"),
                        gp_assoc=detail.get("gp_assoc"),
                        statistical_method=detail.get("statistical_method"),
                        sex=detail.get("sex"),
                        original_id=None,
                        provenance=obj.get("provenance"),
                    )
                )

        return rows

    def protein_mouse_phenotype_converter(self, obj: dict) -> List[mysqlPhenotype]:
        term_name = (obj.get("end_node") or {}).get("name")
        rows: List[mysqlPhenotype] = []
        protein_id = self.resolve_id("protein", obj["start_id"])

        for detail in obj.get("details") or []:
            source = detail.get("source") or DataSourceName.MGI.value
            ptype = "JAX/MGI Human Ortholog Phenotype" if source == DataSourceName.MGI.value else source
            rows.append(
                mysqlPhenotype(
                    ptype=ptype,
                    protein_id=protein_id,
                    nhprotein_id=None,
                    trait=detail.get("trait"),
                    top_level_term_id=detail.get("top_level_term_id"),
                    top_level_term_name=detail.get("top_level_term_name"),
                    term_id=obj["end_id"],
                    term_name=term_name,
                    term_description=None,
                    p_value=detail.get("p_value"),
                    percentage_change=detail.get("percentage_change"),
                    effect_size=detail.get("effect_size"),
                    procedure_name=detail.get("procedure_name"),
                    parameter_name=detail.get("parameter_name"),
                    gp_assoc=detail.get("gp_assoc"),
                    statistical_method=detail.get("statistical_method"),
                    sex=detail.get("sex"),
                    original_id=None,
                    provenance=obj.get("provenance"),
                )
            )

        return rows

    # --- GeneRif ---

    def _iter_publications(self, obj: dict) -> List[dict]:
        return list(obj.get('publications') or [])

    def _group_generifs_from_publications(self, obj: dict) -> List[dict]:
        grouped = {}
        for publication in self._iter_publications(obj):
            gene_id = publication.get('gene_id')
            pmid = publication.get('pmid')
            for gene_rif in publication.get('gene_rifs') or []:
                rif_text = gene_rif.get('text')
                updated_at = gene_rif.get('updated_at')
                if rif_text is None:
                    continue
                key = (gene_id, rif_text)
                entry = grouped.setdefault(
                    key,
                    {
                        'start_id': obj['id'],
                        'gene_id': gene_id,
                        'date': updated_at,
                        'pmids': set(),
                        'end_node': {'text': rif_text},
                        'provenance': obj['provenance'],
                    },
                )
                if updated_at and (entry['date'] is None or updated_at > entry['date']):
                    entry['date'] = updated_at
                if pmid is not None:
                    entry['pmids'].add(int(pmid))

        grouped_rows = []
        for entry in grouped.values():
            entry['pmids'] = sorted(entry['pmids'])
            grouped_rows.append(entry)
        return grouped_rows

    def _generif_hash(self, obj: dict) -> str:
        data_to_hash = f"{obj['start_id']}-{obj['gene_id']}-{obj['date']}-{obj['end_node']['text']}"
        return hashlib.sha256(data_to_hash.encode('utf-8')).hexdigest()

    @staticmethod
    def _coerce_date(value):
        if value is None or value == '':
            return None
        if isinstance(value, date):
            return value
        if isinstance(value, datetime):
            return value.date()
        if isinstance(value, str):
            try:
                return datetime.fromisoformat(value).date()
            except ValueError:
                try:
                    return date.fromisoformat(value)
                except ValueError:
                    return None
        return None

    def generif_converter(self, obj: dict) -> GeneRif:
        return GeneRif(
            id=self._generif_hash(obj),
            protein_id=self.resolve_id('protein', obj['start_id']),
            gene_id=obj['gene_id'],
            date=self._coerce_date(obj['date']),
            text=obj['end_node']['text'],
            provenance=obj['provenance']
        )

    def generif_assoc_converter(self, obj: dict) -> List[GeneRif2Pubmed]:
        hashed_id = self._generif_hash(obj)
        return [
            GeneRif2Pubmed(generif_id=hashed_id, pubmed_id=pmid, provenance=obj['provenance'])
            for pmid in obj['pmids']
        ]

    def p2p_converter(self, obj: dict) -> List[Protein2Pubmed]:
        protein_id = self.resolve_id('protein', obj['start_id'])
        rows = []
        for pmid in obj['pmids']:
            key = (protein_id, str(pmid), obj['gene_id'], 'NCBI')
            if key in self._seen_protein2pubmed:
                continue
            self._seen_protein2pubmed.add(key)
            rows.append(
                Protein2Pubmed(
                    protein_id=protein_id,
                    pubmed_id=pmid,
                    gene_id=obj['gene_id'],
                    source='NCBI',
                    provenance=obj['provenance']
                )
            )
        return rows

    def generif_from_publications_converter(self, obj: dict) -> List[GeneRif]:
        return [self.generif_converter(grouped) for grouped in self._group_generifs_from_publications(obj)]

    def generif_assoc_from_publications_converter(self, obj: dict) -> List[GeneRif2Pubmed]:
        rows: List[GeneRif2Pubmed] = []
        for grouped in self._group_generifs_from_publications(obj):
            rows.extend(self.generif_assoc_converter(grouped))
        return rows

    def p2p_from_publications_converter(self, obj: dict) -> List[Protein2Pubmed]:
        protein_id = self.resolve_id('protein', obj['id'])
        rows = []
        for publication in self._iter_publications(obj):
            pmid = publication.get('pmid')
            source = publication.get('source')
            gene_id = publication.get('gene_id')
            if pmid is None or source is None:
                continue
            key = (protein_id, str(pmid), gene_id, source)
            if key in self._seen_protein2pubmed:
                continue
            self._seen_protein2pubmed.add(key)
            rows.append(
                Protein2Pubmed(
                    protein_id=protein_id,
                    pubmed_id=int(pmid),
                    gene_id=gene_id,
                    source=source,
                    provenance=obj['provenance']
                )
            )
        return rows

    # --- GO ---

    def goterm_converter(self, obj: dict) -> GO:
        def get_namespace(obj: dict):
            if obj['type'] == GoType.Process:
                return 'biological_process'
            if obj['type'] == GoType.Function:
                return 'molecular_function'
            if obj['type'] == GoType.Component:
                return 'cellular_component'
            return 'unknown'
        return GO(
            go_id=obj['id'],
            name=obj['term'],
            namespace=get_namespace(obj),
            def_=obj['definition'],
            provenance=obj['provenance'])

    def goterm_parent_converter(self, obj: dict) -> GOParent:
        return GOParent(
            go_id=obj['start_id'],
            parent_id=obj['end_id'],
            provenance=obj['provenance'])

    def goa_converter(self, obj: dict) -> List[GoA]:
        return [
            GoA(
                protein_id=self.resolve_id('protein', obj['start_id']),
                go_id=obj['end_id'],
                evidence=e['abbreviation'],
                goeco=e['evidence'],
                assigned_by=e['assigned_by'],
                go_term=f"{obj['end_node']['type']}:{obj['end_node']['term']}",
                go_type={'C': 'Component', 'P': 'Process', 'F': 'Function'}.get(obj['end_node']['type'], obj['end_node']['type']),
                go_term_text=obj['end_node']['term'],
                provenance=obj['provenance']
            )
            for e in obj['evidence']]

    # --- Ligand ---

    def ligand_converter(self, obj: dict) -> mysqlLigand:
        xrefs = [EquivalentId.parse(eq_id) for eq_id in (obj.get('xref') or [])]
        ligand_identifier = obj.get('id')
        return mysqlLigand(
            id=self.resolve_id('ligand', ligand_identifier),
            identifier=ligand_identifier,
            name=obj.get('name'),
            isDrug=obj.get('isDrug', False),
            smiles=obj.get('smiles'),
            description=obj.get('description'),
            PubChem=",".join(x.id for x in xrefs if x.type == Prefix.PUBCHEM_COMPOUND),
            ChEMBL=",".join(x.id for x in xrefs if x.type == Prefix.CHEMBL_COMPOUND),
            guide_to_pharmacology=",".join(x.id for x in xrefs if x.type == Prefix.GTOPDB),
            DrugCentral=",".join(x.id for x in xrefs if x.type == Prefix.DrugCentral),
            unii=",".join(x.id for x in xrefs if x.type == Prefix.UNII),
            provenance=obj['provenance']
        )

    def ligand_edge_converter(self, obj: dict) -> List[LigandActivity]:
        activity_objects = []
        ligand_id = self.resolve_id('ligand', obj['end_id'])
        for detail in obj.get('details', []):
            pubmed_ids = list(detail.get('act_pmids') or [])
            if detail.get('moa_pmid'):
                pubmed_ids.append(detail['moa_pmid'])
            activity_objects.append(LigandActivity(
                ncats_ligand_id=ligand_id,
                target_id=self.resolve_id('protein', obj['start_id']),
                act_value=detail.get('act_value'),
                act_type=detail.get('act_type'),
                action_type=detail.get('action_type'),
                reference=detail.get('reference'),
                pubmed_ids="|".join(list(set(str(i) for i in pubmed_ids))),
                provenance=obj['provenance']
            ))
        return activity_objects

    # --- Tissue / Uberon ---

    def uberon_converter(self, obj: dict) -> Uberon:
        return Uberon(
            uid=obj['id'],
            name=obj['name'],
            def_=obj.get('definition'),
            provenance=obj['provenance'],
        )

    def uberon_parent_converter(self, obj: dict) -> UberonParent:
        return UberonParent(
            uid=obj['start_id'],
            parent_id=obj['end_id'],
            provenance=obj['provenance'],
        )

    def tissue_lookup_converter(self, obj: dict) -> List[mysqlTissue]:
        """Populates the tissue name→id lookup table; one row per new unique tissue name."""
        rows = []
        for detail in obj.get('details', []):
            if detail.get('source') == 'GTEx':
                continue
            tissue_name = detail.get('tissue')
            if not tissue_name:
                continue
            is_new = tissue_name not in self.id_mapping.get('tissue', {})
            tissue_id = self.resolve_id('tissue', tissue_name)
            if is_new:
                rows.append(mysqlTissue(id=tissue_id, name=tissue_name))
        return rows

    def expression_converter(self, obj: dict) -> List[Expression]:
        protein_id = self.resolve_id('protein', obj['start_id'])
        resolved_tissue_ontology_id = obj.get('end_id')
        rows = []
        for detail in obj.get('details', []):
            if detail.get('source') == 'GTEx':
                continue
            tissue_name = detail.get('tissue') or ''
            tissue_id = self.resolve_id('tissue', tissue_name) if tissue_name else None
            rows.append(Expression(
                etype=detail['source'],
                protein_id=protein_id,
                source_id=detail.get('source_id') or '',
                tissue=tissue_name,
                tissue_id=tissue_id,
                qual_value=detail.get('qual_value'),
                number_value=detail.get('number_value'),
                expressed=1 if detail.get('expressed') else 0,
                source_rank=detail.get('source_rank'),
                evidence=detail.get('evidence'),
                uberon_id=resolved_tissue_ontology_id or detail.get('source_tissue_id'),
                provenance=obj['provenance'],
            ))
        return rows

    def gtex_converter(self, obj: dict) -> List[Gtex]:
        """One Gtex row per protein-tissue pair, with tpm/tpm_male/tpm_female from sex-split details."""
        protein_id = self.resolve_id('protein', obj['start_id'])
        gtex_details = [d for d in obj.get('details', []) if d.get('source') == 'GTEx']
        if not gtex_details:
            return []

        # Group by source_tissue_id (uberon id of the gtex tissue)
        by_tissue: dict = {}
        for detail in gtex_details:
            tissue_id = detail.get('source_tissue_id') or detail.get('tissue') or ''
            if tissue_id not in by_tissue:
                by_tissue[tissue_id] = {
                    'tissue': detail.get('tissue'),
                    'uberon_id': detail.get('source_tissue_id'),
                    'tpm': None, 'tpm_rank': None,
                    'tpm_male': None, 'tpm_male_rank': None,
                    'tpm_female': None, 'tpm_female_rank': None,
                }
            sex = detail.get('sex')
            if sex == 'male':
                by_tissue[tissue_id]['tpm_male'] = detail.get('number_value')
                by_tissue[tissue_id]['tpm_male_rank'] = detail.get('source_rank')
            elif sex == 'female':
                by_tissue[tissue_id]['tpm_female'] = detail.get('number_value')
                by_tissue[tissue_id]['tpm_female_rank'] = detail.get('source_rank')
            else:
                by_tissue[tissue_id]['tpm'] = detail.get('number_value')
                by_tissue[tissue_id]['tpm_rank'] = detail.get('source_rank')

        return [
            Gtex(
                protein_id=protein_id,
                tissue=d['tissue'],
                tpm=d['tpm'],
                tpm_rank=d['tpm_rank'],
                tpm_male=d['tpm_male'],
                tpm_male_rank=d['tpm_male_rank'],
                tpm_female=d['tpm_female'],
                tpm_female_rank=d['tpm_female_rank'],
                uberon_id=d['uberon_id'],
                provenance=obj['provenance'],
            )
            for d in by_tissue.values()
        ]

    # --- Disease / MONDO ---

    def mondo_table_converter(self, obj: dict) -> Mondo:
        mondoid = obj['id']
        self._known_mondo_ids.add(mondoid)
        return Mondo(
            mondoid=mondoid,
            name=obj.get('name') or mondoid,
            def_=obj.get('mondo_description'),
            comment=obj.get('comment'),
            provenance=obj['provenance'],
        )

    def mondo_xref_converter(self, obj: dict) -> List[MondoXref]:
        rows = []
        seen = set()
        exact_matches = set(obj.get('exact_matches') or [])
        for xref in obj.get('mondo_xrefs') or []:
            if not xref or ':' not in xref or xref in seen:
                continue
            seen.add(xref)
            db, value = xref.split(':', 1)
            rows.append(MondoXref(
                mondoid=obj['id'],
                db=db,
                value=value,
                equiv_to=(xref in exact_matches),
                xref=xref,
                provenance=obj['provenance'],
            ))
        return rows

    def mondo_parent_table_converter(self, obj: dict) -> MondoParent:
        return MondoParent(
            mondoid=obj['start_id'],
            parentid=obj['end_id'],
            provenance=obj['provenance'],
        )

    def do_table_converter(self, obj: dict) -> DO:
        return DO(
            doid=obj['id'],
            name=obj.get('name') or obj['id'],
            def_=obj.get('do_description'),
            provenance=obj['provenance'],
        )

    def do_parent_table_converter(self, obj: dict) -> DOParent:
        return DOParent(
            doid=obj['start_id'],
            parent_id=obj['end_id'],
            provenance=obj['provenance'],
        )

    def mondo_converter(self, obj: dict) -> Optional[Mondo]:
        if not (obj.get('id') or '').startswith('MONDO:'):
            return None
        comments = obj.get('comments') or []
        mondoid = obj['id']
        self._known_mondo_ids.add(mondoid)
        return Mondo(
            mondoid=mondoid,
            name=obj['name'],
            def_=obj.get('mondo_description'),
            comment=comments[0] if comments else None,
            provenance=obj['provenance'],
        )

    def do_converter(self, obj: dict) -> Optional[DO]:
        doid = obj.get('id') or ''
        if not doid.startswith('DOID:'):
            return None
        return DO(
            doid=doid,
            name=obj.get('name') or '',
            def_=obj.get('do_description'),
        )

    def do_parent_converter(self, obj: dict) -> Optional[DOParent]:
        child_id = obj.get('start_id') or ''
        parent_id = obj.get('end_id') or ''
        if not child_id.startswith('DOID:') or not parent_id.startswith('DOID:'):
            return None
        return DOParent(
            doid=child_id,
            parent_id=parent_id,
        )

    def mondo_parent_converter(self, obj: dict) -> MondoParent:
        return MondoParent(
            mondoid=obj['start_id'],
            parentid=obj['end_id'],
            provenance=obj['provenance'],
        )

    def disease_type_converter(self, obj: dict) -> Optional[DiseaseType]:
        """Seeds the disease_type lookup table; one row per new unique source."""
        details = obj.get('details') or []
        source = obj.get('source') or (details[0].get('source') if details else None)
        if not source or source in self._known_disease_types:
            return None
        self._known_disease_types.add(source)
        return DiseaseType(name=source, provenance=obj['provenance'])

    @staticmethod
    def _disease_assoc_key(protein_id: str, disease_id: str, detail: dict, ordinal: int) -> str:
        evidence_terms = "|".join(sorted(detail.get('evidence_terms') or []))
        pmids = "|".join(sorted(str(p) for p in (detail.get('pmids') or [])))
        evidence_codes = "|".join(sorted(detail.get('evidence_codes') or []))
        return "\t".join([
            protein_id or '',
            disease_id or '',
            detail.get('source') or '',
            detail.get('source_id') or '',
            evidence_terms,
            pmids,
            evidence_codes,
            str(ordinal),
        ])

    def _iter_disease_details(self, obj: dict):
        details = obj.get('details') or []
        if details:
            return details
        source = obj.get('source')
        if source:
            return [{
                'source': source,
                'source_id': obj.get('source_id'),
                'evidence_terms': obj.get('evidence_terms') or [],
                'pmids': obj.get('pmids') or [],
                'evidence_codes': obj.get('evidence_codes') or [],
            }]
        return []

    def _disease_name(self, obj: dict) -> str:
        end_id = obj['end_id']
        end_name = (obj.get('end_node') or {}).get('name')
        return end_name or self._disease_name_by_id.get(end_id) or end_id

    def _mondoid_for_disease(self, obj: dict) -> str | None:
        disease_id = obj['id']
        if not disease_id.startswith('IFXDisease:'):
            return disease_id if disease_id.startswith('MONDO:') and disease_id in self._known_mondo_ids else None
        payload = obj.get('disease_harmonizer_payload') or {}
        mondo_ids = {token.strip() for token in payload.get('MONDO_xref', '').split('|') if token.strip()}
        if len(mondo_ids) > 1:
            raise ValueError(f"{disease_id} has multiple exact MONDO aliases: {sorted(mondo_ids)}")
        if not mondo_ids:
            return None
        mondoid = next(iter(mondo_ids))
        if mondoid not in self._known_mondo_ids:
            return None
        return mondoid

    def ncats_disease_converter(self, obj: dict) -> NcatsDisease:
        disease_id = obj['id']
        self._disease_name_by_id[disease_id] = obj.get('name') or disease_id
        mondoid = self._mondoid_for_disease(obj)
        self._disease_mondoid_by_id[disease_id] = mondoid
        payload = obj.get('disease_harmonizer_payload')
        if disease_id.startswith('IFXDisease:') and not isinstance(payload, dict):
            raise ValueError(f"{disease_id} has no DISEASES payload for Pharos conversion")
        rare = obj.get('rare_disease')
        if isinstance(payload, dict):
            raw_rare = payload.get('is_rare')
            if raw_rare not in {'true', 'false'}:
                raise ValueError(f"{disease_id} has invalid DISEASES is_rare {raw_rare!r}")
            rare = raw_rare == 'true'
        return NcatsDisease(
            id=self.resolve_id('ncats_disease', disease_id),
            ifx_id=disease_id if disease_id.startswith('IFXDisease:') else None,
            name=obj.get('name') or disease_id,
            uniprot_description=obj.get('uniprot_description'),
            do_description=obj.get('do_description'),
            mondo_description=obj.get('mondo_description'),
            mondoid=mondoid,
            gard_rare=rare,
            novelty=min(obj.get('novelty') or []) if obj.get('novelty') else None,
            provenance=obj['provenance'],
        )

    def disease_converter(self, obj: dict) -> List[mysqlDisease]:
        resolved_disease_id = obj['end_id']
        if resolved_disease_id and resolved_disease_id.startswith('IFXDisease:'):
            if resolved_disease_id not in self._disease_mondoid_by_id:
                raise ValueError(
                    f"DISEASES node {resolved_disease_id} must be converted before its associations"
                )
            mondoid = self._disease_mondoid_by_id[resolved_disease_id]
        else:
            mondoid = (
                resolved_disease_id
                if resolved_disease_id and resolved_disease_id.startswith('MONDO:') and resolved_disease_id in self._known_mondo_ids
                else None
            )
        disease_name = self._disease_name(obj)
        rows = []
        for ordinal, detail in enumerate(self._iter_disease_details(obj)):
            assoc_key = self._disease_assoc_key(obj['start_id'], resolved_disease_id, detail, ordinal)
            rows.append(mysqlDisease(
                id=self.resolve_id('disease_assoc', assoc_key),
                dtype=detail.get('source') or '',
                protein_id=self.resolve_id('protein', obj['start_id']),
                name=disease_name,
                ncats_name=disease_name,
                did=detail.get('source_id'),
                evidence="|".join(detail.get('evidence_terms') or detail.get('evidence_codes') or []) or None,
                zscore=detail.get('zscore'),
                conf=detail.get('confidence'),
                reference=detail.get('url'),
                drug_name=detail.get('drug_name'),
                source="|".join(detail.get('original_sources') or []) or None,
                mondoid=mondoid,
                provenance=obj['provenance'],
            ))
        return rows

    def ncats_d2da_converter(self, obj: dict) -> List[NcatsD2DA]:
        resolved_disease_id = obj['end_id']
        ncats_disease_id = self.resolve_id('ncats_disease', resolved_disease_id)
        links = []
        for ordinal, detail in enumerate(self._iter_disease_details(obj)):
            assoc_key = self._disease_assoc_key(obj['start_id'], resolved_disease_id, detail, ordinal)
            links.append(NcatsD2DA(
                ncats_disease_id=ncats_disease_id,
                disease_assoc_id=self.resolve_id('disease_assoc', assoc_key),
                direct=1,
                provenance=obj['provenance'],
            ))
        return links

    def tinx_importance_converter(self, obj: dict) -> List[TinxImportance]:
        best_doid = None
        best_score = None
        for detail in self._iter_disease_details(obj):
            importance = detail.get('importance') or []
            doid = detail.get('doid') or detail.get('source_id')
            if not importance:
                continue
            score = max(importance)
            if best_score is None or score > best_score:
                best_score = score
                best_doid = doid
        if best_score is None:
            return []
        return [TinxImportance(
            ncats_disease_id=self.resolve_id('ncats_disease', obj['end_id']),
            protein_id=self.resolve_id('protein', obj['start_id']),
            doid=best_doid,
            score=best_score,
        )]

    # --- TIGA ---

    def tiga_converter(self, obj: dict) -> List[mysqlTiga]:
        trait = obj.get('end_node') or {}
        if not trait:
            return []
        rows = []
        disease_ids = obj.get('disease_ids') or [None]
        for detail in obj.get('details') or []:
            for disease_id in disease_ids:
                rows.append(mysqlTiga(
                    protein_id=self.resolve_id('protein', obj['start_id']),
                    ensg=detail.get('ensg') or '',
                    efoid=obj.get('end_id') or '',
                    trait=trait.get('name') or '',
                    n_study=detail.get('n_study'),
                    n_snp=detail.get('n_snp'),
                    n_snpw=detail.get('n_snpw'),
                    geneNtrait=detail.get('geneNtrait'),
                    geneNstudy=detail.get('geneNstudy'),
                    traitNgene=detail.get('traitNgene'),
                    traitNstudy=detail.get('traitNstudy'),
                    pvalue_mlog_median=detail.get('pvalue_mlog_median'),
                    pvalue_mlog_max=detail.get('pvalue_mlog_max'),
                    or_median=detail.get('or_median'),
                    n_beta=detail.get('n_beta'),
                    study_N_mean=int(detail['study_N_mean']) if detail.get('study_N_mean') is not None else None,
                    rcras=detail.get('rcras'),
                    meanRank=detail.get('meanRank'),
                    meanRankScore=detail.get('meanRankScore'),
                    ncats_disease_id=self.resolve_id('ncats_disease', disease_id) if disease_id else None,
                    provenance=obj['provenance'],
                ))
        return rows

    def tiga_provenance_converter(self, obj: dict) -> List[TigaProvenance]:
        if not obj.get('end_node'):
            return []
        rows = []
        for edge_detail in obj.get('details') or []:
            for prov_detail in edge_detail.get('provenance_details') or []:
                study_acc = prov_detail.get('study_acc')
                pubmedid = prov_detail.get('pubmedid')
                if not study_acc or pubmedid is None:
                    continue
                rows.append(TigaProvenance(
                    ensg=edge_detail.get('ensg') or '',
                    efoid=obj.get('end_id') or '',
                    study_acc=study_acc,
                    pubmedid=int(pubmedid),
                    provenance=obj['provenance'],
                ))
        return rows

    # --- Pathway ---

    def pathway_converter(self, obj: dict) -> mysqlPathway:
        protein_id = self.resolve_id('protein', obj['start_id'])
        end = obj['end_node']
        return mysqlPathway(
            target_id=protein_id,
            protein_id=protein_id,
            pwtype=end.get('type'),
            id_in_source=end.get('source_id'),
            name=end.get('name') or '',
            url=end.get('url'),
            provenance=obj['provenance'],
        )

    # --- PPI ---

    @staticmethod
    def _ppi_source_label(source: str) -> str:
        if not source:
            return source
        label = source.split('\t', 1)[0]
        return {'STRING': 'StringDB'}.get(label, label)

    @staticmethod
    def _max_or_value(value):
        if isinstance(value, list):
            return max(value) if value else None
        return value

    @staticmethod
    def _join_list(value):
        if isinstance(value, list):
            normalized = [str(v) for v in value if v is not None and str(v) != ""]
            return "|".join(dict.fromkeys(normalized)) if normalized else None
        return value

    def ppi_converter(self, obj: dict) -> List[mysqlPPI]:
        protein_id = self.resolve_id('protein', obj['start_id'])
        other_id = self.resolve_id('protein', obj['end_id'])

        source_labels = sorted({
            self._ppi_source_label(source)
            for source in (obj.get('sources') or [])
            if source
        })

        shared = dict(
            ppitypes=",".join(source_labels) if source_labels else self._ppi_source_label(obj['provenance']),
            p_int=self._max_or_value(obj.get('p_int')),
            p_ni=self._max_or_value(obj.get('p_ni')),
            p_wrong=self._max_or_value(obj.get('p_wrong')),
            evidence=self._join_list(obj.get('pmids')) or obj.get('evidence'),
            interaction_type=self._join_list(obj.get('interaction_type')),
            score=self._max_or_value(obj.get('score')),
        )

        return [
            mysqlPPI(
                protein_id=protein_id,
                other_id=other_id,
                **shared,
            ),
            mysqlPPI(
                protein_id=other_id,
                other_id=protein_id,
                **shared,
            ),
        ]

    # --- Viral PPI ---

    @staticmethod
    def _viral_protein_pk(graph_id: str) -> int:
        return int(graph_id.rsplit(":", 1)[1])

    @staticmethod
    def _virus_taxid(graph_id: str) -> str:
        return graph_id.split(":", 1)[1]

    def virus_converter(self, obj: dict) -> mysqlVirus:
        return mysqlVirus(
            virusTaxid=self._virus_taxid(obj["id"]),
            nucleic1=obj.get("nucleic1"),
            nucleic2=obj.get("nucleic2"),
            order=obj.get("order"),
            family=obj.get("family"),
            subfamily=obj.get("subfamily"),
            genus=obj.get("genus"),
            species=obj.get("species"),
            name=obj.get("name"),
        )

    def viral_protein_converter(self, obj: dict) -> mysqlViralProtein:
        viral_protein_to_virus = self.id_mapping.get("viral_protein_to_virus", {})
        return mysqlViralProtein(
            id=self._viral_protein_pk(obj["id"]),
            name=obj.get("name"),
            ncbi=obj.get("ncbi"),
            virus_id=viral_protein_to_virus.get(obj["id"]),
        )

    @staticmethod
    def _viral_ppi_details(obj: dict) -> List[ViralPPIDetail]:
        details = obj.get("details") or []
        if details:
            return [ViralPPIDetail.from_dict(detail) for detail in details]

        # Backward-compatible fallback for older graph builds before PHIPSTER details
        legacy_source = obj.get("data_source")
        legacy_final_lr = obj.get("final_lr")
        legacy_pdb_ids = obj.get("pdb_ids") or []
        legacy_high_confidence = obj.get("high_confidence")
        if (
            legacy_source is None
            and legacy_final_lr is None
            and not legacy_pdb_ids
            and legacy_high_confidence is None
        ):
            return []
        return [
            ViralPPIDetail(
                source=legacy_source,
                source_protein_id=obj.get("start_id"),
                final_lr=legacy_final_lr,
                pdb_ids=legacy_pdb_ids,
                high_confidence=legacy_high_confidence,
            )
        ]

    def viral_ppi_converter(self, obj: dict) -> List[mysqlViralPPI]:
        protein_id = self.resolve_id("protein", obj["start_id"])
        viral_protein_id = self._viral_protein_pk(obj["end_id"])
        rows = []
        for detail in self._viral_ppi_details(obj):
            rows.append(
                mysqlViralPPI(
                    viral_protein_id=viral_protein_id,
                    protein_id=protein_id,
                    dataSource=detail.source or self._ppi_source_label(obj.get("provenance")),
                    finalLR=detail.final_lr,
                    pdbIDs=self._join_list(detail.pdb_ids),
                    highConfidence=detail.high_confidence,
                )
            )
        return rows

    # --- Panther ---

    def panther_class_converter(self, obj: dict) -> mysqlPantherClass:
        return mysqlPantherClass(
            id=self.resolve_id('panther_class', obj['id']),
            pcid=obj.get('source_id') or obj['id'],
            parent_pcids=obj.get('parent_pcids'),
            name=obj.get('name') or '',
            description=obj.get('description'),
        )

    def p2pc_converter(self, obj: dict) -> P2PC:
        return P2PC(
            panther_class_id=self.resolve_id('panther_class', obj['end_id']),
            protein_id=self.resolve_id('protein', obj['start_id']),
        )

    # --- DTO ---

    def dto_converter(self, obj: dict) -> mysqlDTO:
        return mysqlDTO(
            dtoid=obj['id'],
            name=obj.get('name') or '',
            parent_id=None,
            def_=obj.get('description'),
        )

    def dto_parent_converter(self, obj: dict) -> DTOParent:
        return DTOParent(
            dtoid=obj['start_id'],
            parent_id=obj['end_id'],
        )

    def p2dto_converter(self, obj: dict) -> P2DTO:
        return P2DTO(
            dtoid=obj['end_id'],
            protein_id=self.resolve_id('protein', obj['start_id']),
        )

    # --- Keyword ---

    def keyword_xref_converter(self, obj: dict) -> Xref:
        end = obj['end_node']
        return Xref(
            protein_id=self.resolve_id('protein', obj['start_id']),
            xtype='UniProt Keyword',
            value=end.get('source_id') or end.get('value') or end.get('name') or '',
            xtra=end.get('value') or end.get('name'),
            provenance=obj['provenance'],
        )

    # --- Simple Linkouts ---

    def affiliate_converter(self, obj: dict) -> Affiliate:
        source = obj.get("source")
        return Affiliate(
            source=source,
            display_name=obj.get("name") or source,
            description=obj.get("description") or "",
            link_count=obj.get("link_count") or 0,
        )

    @staticmethod
    def _pick_extlink_detail(obj: dict) -> dict:
        details = obj.get("details") or []
        if details:
            return sorted(details, key=lambda detail: detail.get("source_id") or "")[0]
        return {
            "url": obj.get("url"),
            "source_id": obj.get("source_id"),
            "source_id_type": obj.get("source_id_type"),
        }

    def extlink_converter(self, obj: dict) -> Optional[ExtLink]:
        detail = self._pick_extlink_detail(obj)
        url = detail.get("url")
        if not url:
            return None
        return ExtLink(
            protein_id=self.resolve_id("protein", obj["start_id"]),
            source=obj.get("source"),
            url=url,
        )
