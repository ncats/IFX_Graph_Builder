from pathlib import Path
import csv
import gzip
from typing import Dict, Generator, Iterable, List, Optional, Set, Tuple, Union

import ijson
from rdflib import Graph, URIRef

from src.constants import DataSourceName, Prefix
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite_harmonization import (
    MetaboliteIdentifier,
    ProteinIdentifier,
    RheaMetaboliteReactionEdge,
    RheaProteinReactionEdge,
    RheaReaction,
    RheaReactionDirectionEdge,
    RheaReactionReactionClassEdge,
    RheaReactionClass,
)


RHEA = "http://rdf.rhea-db.org/"
RDFS = "http://www.w3.org/2000/01/rdf-schema#"


class RheaReactionAdapter(InputAdapter):
    def __init__(
        self,
        data_source=None,
        uniprot_data_source=None,
        rdf_file: Optional[str] = None,
        directions_file: Optional[str] = None,
        rhea2ec_file: Optional[str] = None,
        rhea2uniprot_sprot_file: Optional[str] = None,
        rhea2uniprot_trembl_file: Optional[str] = None,
        uniprot_human_file: Optional[str] = None,
        active_only: bool = True,
        filter_human_proteins: bool = True,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            try:
                rdf_file = str(data_source.file("rhea.rdf.gz"))
            except FileNotFoundError:
                rdf_file = str(data_source.file("rhea.rdf"))
            directions_file = str(data_source.file("rhea-directions.tsv"))
            rhea2ec_file = str(data_source.file("rhea2ec.tsv"))
            rhea2uniprot_sprot_file = str(data_source.file("rhea2uniprot_sprot.tsv"))
            rhea2uniprot_trembl_file = str(data_source.file("rhea2uniprot_trembl.tsv.gz"))
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(version=None, version_date=None, download_date=None)

        if uniprot_data_source is not None:
            uniprot_human_file = str(uniprot_data_source.file("uniprot-human.json.gz"))

        required = {
            "rdf_file": rdf_file,
            "directions_file": directions_file,
            "rhea2ec_file": rhea2ec_file,
            "rhea2uniprot_sprot_file": rhea2uniprot_sprot_file,
            "rhea2uniprot_trembl_file": rhea2uniprot_trembl_file,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise ValueError("RheaReactionAdapter missing required files: " + ", ".join(missing))
        if filter_human_proteins and uniprot_human_file is None:
            raise ValueError("RheaReactionAdapter requires uniprot_data_source or uniprot_human_file for human filtering")

        self.rdf_file = Path(rdf_file)
        self.directions_file = Path(directions_file)
        self.rhea2ec_file = Path(rhea2ec_file)
        self.rhea2uniprot_sprot_file = Path(rhea2uniprot_sprot_file)
        self.rhea2uniprot_trembl_file = Path(rhea2uniprot_trembl_file)
        self.uniprot_human_file = None if uniprot_human_file is None else Path(uniprot_human_file)
        self.active_only = active_only
        self.filter_human_proteins = filter_human_proteins
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.RHEA

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(
        self,
    ) -> Generator[
        List[
            Union[
                RheaReaction,
                MetaboliteIdentifier,
                ProteinIdentifier,
                RheaReactionDirectionEdge,
                RheaMetaboliteReactionEdge,
                RheaProteinReactionEdge,
                RheaReactionReactionClassEdge,
            ]
        ],
        None,
        None,
    ]:
        reactions, master_to_variants = self._load_reactions()
        yield from self._yield_batches(self._iter_reaction_nodes(reactions))
        yield from self._yield_batches(self._iter_metabolite_nodes(reactions))
        yield from self._yield_batches(self._iter_protein_nodes(reactions))
        yield from self._yield_batches(self._iter_direction_edges(master_to_variants, reactions))
        yield from self._yield_batches(self._iter_metabolite_reaction_edges(reactions))
        yield from self._yield_batches(self._iter_protein_reaction_edges(reactions))
        yield from self._yield_batches(self._iter_reaction_class_edges(reactions))

    def _yield_batches(self, objects):
        batch = []
        for obj in objects:
            batch.append(obj)
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    def _load_reactions(self) -> Tuple[Dict[str, Dict], Dict[str, Dict[str, str]]]:
        graph = Graph()
        if self.rdf_file.suffix == ".gz":
            with gzip.open(self.rdf_file, "rb") as handle:
                graph.parse(handle, format="application/rdf+xml")
        else:
            graph.parse(self.rdf_file, format="application/rdf+xml")
        reactions = self._parse_reaction_rdf(graph)
        master_to_variants = self._parse_direction_rows()
        self._apply_direction_rows(reactions, master_to_variants)
        self._apply_ec_rows(reactions)

        if self.active_only:
            reactions = {
                reaction_id: reaction
                for reaction_id, reaction in reactions.items()
                if reaction.get("status") == "Approved"
            }
        if self.max_records is not None:
            keep = set(sorted(reactions, key=_rhea_sort_key)[: self.max_records])
            reactions = {reaction_id: reaction for reaction_id, reaction in reactions.items() if reaction_id in keep}
            master_to_variants = {
                master_id: variants
                for master_id, variants in master_to_variants.items()
                if master_id in keep or any(variant_id in keep for variant_id in variants.values())
            }
        self._apply_protein_rows(reactions)
        return reactions, master_to_variants

    def _parse_reaction_rdf(self, graph: Graph) -> Dict[str, Dict]:
        reaction_subjects = set()
        subclass_predicate = URIRef(RDFS + "subClassOf")
        for class_name in ["Reaction", "DirectionalReaction", "BidirectionalReaction"]:
            reaction_subjects.update(graph.subjects(predicate=subclass_predicate, object=URIRef(RHEA + class_name)))

        reactions = {}
        for subject in reaction_subjects:
            reaction_id = _first_text(graph.objects(subject, URIRef(RHEA + "accession")))
            if reaction_id is None:
                continue
            reaction_id = _normalize_rhea_id(reaction_id)
            reactions[reaction_id] = {
                "id": reaction_id,
                "source_uri": str(subject),
                "source_id": reaction_id,
                "label": _first_text(graph.objects(subject, URIRef(RDFS + "label"))),
                "equation": _first_text(graph.objects(subject, URIRef(RHEA + "equation"))),
                "html_equation": _first_text(graph.objects(subject, URIRef(RHEA + "htmlEquation"))),
                "status": _status_value(_first_text(graph.objects(subject, URIRef(RHEA + "status")))),
                "is_transport": _first_text(graph.objects(subject, URIRef(RHEA + "isTransport"))) == "true",
                "left": self._participants_for_side(graph, subject, "L"),
                "right": self._participants_for_side(graph, subject, "R"),
                "ec_ids": set(),
                "proteins": {},
                "master_id": None,
                "direction": None,
            }
        return reactions

    def _participants_for_side(self, graph: Graph, reaction_subject, side_suffix: str) -> List[Dict]:
        side_subject = URIRef(str(reaction_subject) + "_" + side_suffix)
        participants_by_uri = {}
        for predicate, participant in graph.predicate_objects(side_subject):
            predicate_text = str(predicate)
            if not predicate_text.startswith(RHEA + "contains"):
                continue
            participant_uri = str(participant)
            coefficient = _contains_coefficient(predicate_text)
            existing = participants_by_uri.get(participant_uri)
            if existing is not None and existing.get("coefficient") is not None:
                continue
            participants_by_uri[participant_uri] = self._participant_record(graph, participant, coefficient)
        return list(participants_by_uri.values())

    def _participant_record(self, graph: Graph, participant, coefficient: Optional[int]) -> Dict:
        compound = next(graph.objects(participant, URIRef(RHEA + "compound")), None)
        if compound is None:
            return {"participant_uri": str(participant), "coefficient": coefficient}
        accession = _first_text(graph.objects(compound, URIRef(RHEA + "accession")))
        return {
            "participant_uri": str(participant),
            "compound_uri": str(compound),
            "coefficient": coefficient,
            "rhea_accession": accession,
            "metabolite_id": _metabolite_id(accession),
            "name": _first_text(graph.objects(compound, URIRef(RHEA + "name"))),
            "html_name": _first_text(graph.objects(compound, URIRef(RHEA + "htmlName"))),
            "formula": _first_text(graph.objects(compound, URIRef(RHEA + "formula"))),
        }

    def _parse_direction_rows(self) -> Dict[str, Dict[str, str]]:
        rows = {}
        with self.directions_file.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                master_id = _normalize_rhea_id(row["RHEA_ID_MASTER"])
                rows[master_id] = {
                    "UN": _normalize_rhea_id(row["RHEA_ID_MASTER"]),
                    "LR": _normalize_rhea_id(row["RHEA_ID_LR"]),
                    "RL": _normalize_rhea_id(row["RHEA_ID_RL"]),
                    "BD": _normalize_rhea_id(row["RHEA_ID_BI"]),
                }
        return rows

    @staticmethod
    def _apply_direction_rows(reactions: Dict[str, Dict], master_to_variants: Dict[str, Dict[str, str]]) -> None:
        for master_id, variants in master_to_variants.items():
            master = reactions.get(master_id)
            if master is None:
                continue
            for direction, reaction_id in variants.items():
                reaction = reactions.get(reaction_id)
                if reaction is None:
                    continue
                reaction["master_id"] = master_id
                reaction["direction"] = direction
                if direction == "RL":
                    reaction["left"] = list(master["right"])
                    reaction["right"] = list(master["left"])
                elif direction != "UN":
                    reaction["left"] = list(master["left"])
                    reaction["right"] = list(master["right"])
            master["master_id"] = master_id
            master["direction"] = "UN"

    def _apply_ec_rows(self, reactions: Dict[str, Dict]) -> None:
        with self.rhea2ec_file.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                reaction_id = _normalize_rhea_id(row["RHEA_ID"])
                reaction = reactions.get(reaction_id)
                if reaction is not None:
                    reaction["ec_ids"].add(row["ID"].strip())

    def _apply_protein_rows(self, reactions: Dict[str, Dict]) -> None:
        human_uniprot = self._human_uniprot_accessions() if self.filter_human_proteins else None
        for path, is_reviewed in [
            (self.rhea2uniprot_sprot_file, True),
            (self.rhea2uniprot_trembl_file, False),
        ]:
            with _open_text(path) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    reaction_id = _normalize_rhea_id(row["RHEA_ID"])
                    reaction = reactions.get(reaction_id)
                    if reaction is None:
                        continue
                    accession = (row.get("ID") or "").strip()
                    if not accession:
                        continue
                    if human_uniprot is not None and accession not in human_uniprot:
                        continue
                    # Membership filtering is not identity resolution. Retain the
                    # accession reported by Rhea, including secondary accessions.
                    protein_id = f"{Prefix.UniProtKB}:{accession}"
                    existing = reaction["proteins"].get(protein_id)
                    reaction["proteins"][protein_id] = {
                        "id": protein_id,
                        "is_reviewed": is_reviewed if existing is None else existing["is_reviewed"] or is_reviewed,
                        "source_file": path.name,
                    }

    def _human_uniprot_accessions(self) -> set[str]:
        if self.uniprot_human_file is None:
            return set()
        accessions = set()
        with gzip.open(self.uniprot_human_file, "rb") as handle:
            for record in ijson.items(handle, "results.item"):
                primary = record.get("primaryAccession")
                if not primary:
                    continue
                accessions.add(primary)
                for secondary in record.get("secondaryAccessions") or []:
                    accessions.add(secondary)
        return accessions

    @staticmethod
    def _iter_reaction_nodes(reactions: Dict[str, Dict]) -> Iterable[RheaReaction]:
        for reaction_id in sorted(reactions, key=_rhea_sort_key):
            reaction = reactions[reaction_id]
            yield RheaReaction(
                id=reaction_id,
                source_id=reaction["source_id"],
                master_id=reaction.get("master_id"),
                direction=reaction.get("direction"),
                status=reaction.get("status"),
                is_transport=reaction.get("is_transport"),
                label=reaction.get("label"),
                equation=reaction.get("equation"),
                html_equation=reaction.get("html_equation"),
            )

    @staticmethod
    def _iter_metabolite_nodes(reactions: Dict[str, Dict]) -> Iterable[MetaboliteIdentifier]:
        emitted = set()
        for reaction in reactions.values():
            for participant in [*reaction.get("left", []), *reaction.get("right", [])]:
                metabolite_id = participant.get("metabolite_id")
                if metabolite_id is None or metabolite_id in emitted:
                    continue
                emitted.add(metabolite_id)
                yield MetaboliteIdentifier(id=metabolite_id)

    @staticmethod
    def _iter_protein_nodes(reactions: Dict[str, Dict]) -> Iterable[ProteinIdentifier]:
        emitted = {}
        for reaction in reactions.values():
            for protein_id, protein in reaction["proteins"].items():
                emitted[protein_id] = emitted.get(protein_id, False) or protein["is_reviewed"]
        for protein_id in sorted(emitted):
            yield ProteinIdentifier(id=protein_id, is_reviewed=emitted[protein_id])

    @staticmethod
    def _iter_direction_edges(
        master_to_variants: Dict[str, Dict[str, str]],
        reactions: Dict[str, Dict],
    ) -> Iterable[RheaReactionDirectionEdge]:
        for master_id in sorted(master_to_variants, key=_rhea_sort_key):
            if master_id not in reactions:
                continue
            for direction in ["LR", "RL", "BD"]:
                reaction_id = master_to_variants[master_id].get(direction)
                if reaction_id not in reactions:
                    continue
                yield RheaReactionDirectionEdge(
                    start_node=RheaReaction(id=master_id),
                    end_node=RheaReaction(id=reaction_id),
                    source_field="rhea-directions.tsv",
                    variant_direction=direction,
                )

    @staticmethod
    def _iter_metabolite_reaction_edges(reactions: Dict[str, Dict]) -> Iterable[RheaMetaboliteReactionEdge]:
        for reaction_id in sorted(reactions, key=_rhea_sort_key):
            reaction = reactions[reaction_id]
            for side in ["left", "right"]:
                for participant in reaction.get(side, []):
                    metabolite_id = participant.get("metabolite_id")
                    if metabolite_id is None:
                        continue
                    yield RheaMetaboliteReactionEdge(
                        start_node=MetaboliteIdentifier(id=metabolite_id),
                        end_node=RheaReaction(id=reaction_id),
                        source_id=metabolite_id,
                        side=side,
                        coefficient=participant.get("coefficient"),
                        participant_uri=participant.get("participant_uri"),
                        compound_uri=participant.get("compound_uri"),
                        rhea_accession=participant.get("rhea_accession"),
                        name=participant.get("name"),
                        html_name=participant.get("html_name"),
                        formula=participant.get("formula"),
                        is_cofactor=None,
                    )

    @staticmethod
    def _iter_protein_reaction_edges(reactions: Dict[str, Dict]) -> Iterable[RheaProteinReactionEdge]:
        for reaction_id in sorted(reactions, key=_rhea_sort_key):
            reaction = reactions[reaction_id]
            for protein_id, protein in sorted(reaction["proteins"].items()):
                yield RheaProteinReactionEdge(
                    start_node=ProteinIdentifier(id=protein_id, is_reviewed=protein["is_reviewed"]),
                    end_node=RheaReaction(id=reaction_id),
                    source_field="rhea2uniprot",
                    source_id=protein_id,
                    source_file=protein.get("source_file"),
                )

    @staticmethod
    def _iter_reaction_class_edges(reactions: Dict[str, Dict]) -> Iterable[RheaReactionReactionClassEdge]:
        for reaction_id in sorted(reactions, key=_rhea_sort_key):
            for ec_id in sorted(reactions[reaction_id]["ec_ids"], key=_ec_sort_key):
                yield RheaReactionReactionClassEdge(
                    start_node=RheaReaction(id=reaction_id),
                    end_node=RheaReactionClass(id=f"EC:{ec_id}"),
                    source_field="rhea2ec.tsv",
                )


def _first_text(values) -> Optional[str]:
    for value in values:
        return str(value).strip()
    return None


def _status_value(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    return value.rsplit("/", 1)[-1]


def _normalize_rhea_id(value: str) -> str:
    text = str(value).strip()
    if ":" in text:
        prefix, suffix = text.split(":", 1)
        if prefix.upper() == "RHEA":
            return f"RHEA:{suffix}"
    return f"RHEA:{text}"


def _metabolite_id(accession: Optional[str]) -> Optional[str]:
    if accession is None:
        return None
    prefix, _, suffix = accession.strip().partition(":")
    if not suffix:
        return None
    prefix_upper = prefix.upper()
    if prefix_upper == "CHEBI":
        return f"CHEBI:{suffix}"
    if prefix_upper == "GENERIC":
        return f"RHEA.COMP:{suffix}"
    if prefix_upper == "POLYMER":
        return f"RHEA.POLYMER:{suffix}"
    return f"{prefix_upper}:{suffix}"


def _contains_coefficient(predicate: str) -> Optional[int]:
    suffix = predicate.rsplit("contains", 1)[-1]
    if suffix.isdigit():
        return int(suffix)
    return None


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _rhea_sort_key(rhea_id: str):
    suffix = rhea_id.split(":", 1)[1]
    return int(suffix) if suffix.isdigit() else suffix


def _ec_sort_key(ec_id: str):
    key = []
    for part in ec_id.split("."):
        if part == "-":
            key.append((-1, 0))
        elif part.startswith("n") and part[1:].isdigit():
            key.append((1, int(part[1:])))
        else:
            key.append((0, int(part)))
    return key
