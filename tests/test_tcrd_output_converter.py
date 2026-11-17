from src.interfaces.input_adapter import InputAdapter
from src.models.disease import Disease, DiseaseAssociationDetail, GeneDiseaseEdge, ProteinDiseaseEdge
from src.models.gene import Gene
from src.models.keyword import ProteinKeywordEdge
from src.models.ligand import Ligand
from src.models.protein import Protein
from src.output_adapters.sql_converters.tcrd import TCRDOutputConverter
from src.shared.sqlalchemy_tables.pharos_tables_new import DOParent, EntityXref, Xref


class _FakeQuery:
    def __init__(self, rows):
        self._rows = rows

    def filter(self, *_args, **_kwargs):
        return self

    def join(self, *_args, **_kwargs):
        return self

    def all(self):
        return self._rows


class _FakeSession:
    def query(self, *columns):
        key = tuple(getattr(column, "name", str(column)) for column in columns)
        if key == ("id", "ifx_id"):
            return _FakeQuery([
                (123, "IFX123"),
            ])
        if key == ("id", "uniprot"):
            return _FakeQuery([
                (123, "Q9Y345"),
            ])
        if key == ("id", "identifier"):
            return _FakeQuery([
                (77, "LIGAND:EXISTING"),
            ])
        if key == ("id", "name"):
            return _FakeQuery([])
        if key == ("id", "pcid"):
            return _FakeQuery([])
        if key == ("name",):
            return _FakeQuery([])
        if key == ("id", "mondoid"):
            return _FakeQuery([])
        if key == ("mondoid",):
            return _FakeQuery([])
        if key == ("geneid", "id"):
            return _FakeQuery([])
        if key == ("sym", "id"):
            return _FakeQuery([])
        raise AssertionError(f"Unexpected query columns: {key}")


def test_tcrd_output_converter_preloads_protein_ids():
    converter = TCRDOutputConverter()

    converter.preload_id_mappings(_FakeSession())

    assert converter.id_mapping["protein"] == {"IFX123": 123}
    assert converter.id_mapping["target_by_uniprot"] == {"Q9Y345": 123}
    assert converter.id_mapping["ligand"] == {"LIGAND:EXISTING": 77}


def test_drgc_resource_converter_uses_target_mapping_by_uniprot():
    converter = TCRDOutputConverter()
    converter._target_id_by_uniprot = {"Q9Y345": 123}

    row = converter.drgc_resource_converter({
        "id": "https://repo.metadatacenter.org/template-instances/example",
        "uniprot_id": "Q9Y345",
        "rssid": "https://repo.metadatacenter.org/template-instances/example",
        "resource_type": "Genetic Construct",
        "json": "{\"resource\": []}",
        "legacy_target_id": 999,
        "provenance": "Pharos 3.19\t3.19\t2024-02-15\tNone",
    })

    assert row.target_id == 123
    assert row.rssid == "https://repo.metadatacenter.org/template-instances/example"
    assert row.resource_type == "Genetic Construct"


def test_tdl_info_converter_emits_surechembl_patent_family_total():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.tdl_info_converter({
        "id": "IFX123",
        "patent_family_mentions": [
            "2020:1001",
            "2020:1001",
            "2021:1002",
        ],
    })

    patent_rows = [row for row in rows if row.itype == "SureChEMBL Patent Family Count"]

    assert len(patent_rows) == 1
    assert patent_rows[0].protein_id == 123
    assert patent_rows[0].integer_value == 2


def test_tdl_info_converter_emits_pubtator_score():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.tdl_info_converter({
        "id": "IFX123",
        "pt_score": [4.25],
    })

    by_type = {row.itype: row for row in rows}
    assert by_type["PubTator Score"].protein_id == 123
    assert by_type["PubTator Score"].number_value == 4.25


def test_patent_count_converter_groups_unique_families_by_year():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.patent_count_converter({
        "id": "IFX123",
        "patent_family_mentions": [
            "2020:1001",
            "2020:1001",
            "2020:1002",
            "2021:1002",
        ],
    })

    assert [(row.protein_id, row.year, row.count) for row in rows] == [
        (123, 2020, 2),
        (123, 2021, 1),
    ]


def test_ptscore_converter_emits_yearly_rows():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.ptscore_converter({
        "id": "IFX123",
        "pt_score_by_year": [
            {"year": 2020, "score": 1.5},
            {"year": 2021, "score": 2.0},
        ],
    })

    assert [(row.protein_id, row.year, float(row.score)) for row in rows] == [
        (123, 2020, 1.5),
        (123, 2021, 2.0),
    ]


def test_drgc_resource_converter_skips_missing_destination_target():
    converter = TCRDOutputConverter()

    row = converter.drgc_resource_converter({
        "id": "rssid",
        "uniprot_id": "Q5JXX5",
        "rssid": "rssid",
        "resource_type": "Mouse",
        "json": "{}",
        "legacy_target_id": 1234,
        "provenance": "Pharos 3.19\t3.19\t2024-02-15\tNone",
    })

    assert row is None


def test_nhprotein_converter_projects_ortholog_gene_into_gene_centric_row():
    converter = TCRDOutputConverter()

    row = converter.nhprotein_converter({
        "id": "NCBIGene:13489",
        "species": "10090",
        "symbol": "Drd2",
        "name": "dopamine receptor D2",
        "entrez_gene_id": "13489",
        "xref": ["NCBIGene:13489", "MGI:94924", "ENSEMBL:ENSMUSG00000032259"],
        "provenance": "HCOP\tlatest\t2026-05-04\t2026-05-05",
    })

    assert row.id == 1
    assert row.uniprot is None
    assert row.name == "Drd2"
    assert row.description == "dopamine receptor D2"
    assert row.sym == "Drd2"
    assert row.species == "Mus musculus"
    assert row.taxid == 10090
    assert row.geneid == 13489


def test_phenotype_converter_emits_one_row_per_detail_per_linked_protein():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {
        "IFX1": 101,
        "IFX2": 202,
    }

    rows = converter.phenotype_converter({
        "start_id": "NCBIGene:13489",
        "end_id": "MP:0005635",
        "protein_links": [
            {"protein_node_id": "IFX1", "support_sources": ["OMA", "EggNOG"]},
            {"protein_node_id": "IFX2", "support_sources": ["Inparanoid", "OMA", "HGNC"]},
        ],
        "end_node": {"name": "decreased circulating bilirubin level"},
        "details": [
            {
                "source": "IMPC",
                "top_level_term_id": "MP:0005376",
                "top_level_term_name": "homeostasis/metabolism phenotype",
                "p_value": 3.97e-09,
                "percentage_change": "-143.75",
                "effect_size": "1.54",
                "procedure_name": "Clinical Chemistry",
                "parameter_name": "Total bilirubin",
                "gp_assoc": True,
                "statistical_method": "GLS",
                "sex": "male",
            },
            {
                "source": "IMPC",
                "top_level_term_id": "MP:0005376",
                "top_level_term_name": "homeostasis/metabolism phenotype",
                "p_value": 1.0e-06,
                "percentage_change": "-120",
                "effect_size": "1.2",
                "procedure_name": "Clinical Chemistry",
                "parameter_name": "Total bilirubin",
                "gp_assoc": True,
                "statistical_method": "GLS",
                "sex": "female",
            },
        ],
        "provenance": "IMPC\tgenotype-phenotype-assertions-IMPC.csv.gz\t2026-03-16\t2026-05-05",
    })

    assert len(rows) == 4
    assert {row.protein_id for row in rows} == {101, 202}
    assert {row.sex for row in rows} == {"male", "female"}
    assert {row.nhprotein_id for row in rows} == {1}
    assert all(row.ptype == "IMPC" for row in rows)
    assert all(row.term_id == "MP:0005635" for row in rows)
    assert all(row.term_name == "decreased circulating bilirubin level" for row in rows)


def test_phenotype_converter_keeps_nhprotein_only_when_no_human_links_present():
    converter = TCRDOutputConverter()

    rows = converter.phenotype_converter({
        "start_id": "NCBIGene:13489",
        "end_id": "MP:0005635",
        "protein_links": [],
        "end_node": {"name": "decreased circulating bilirubin level"},
        "details": [{
            "source": "IMPC",
            "procedure_name": "Clinical Chemistry",
            "parameter_name": "Total bilirubin",
        }],
        "provenance": "IMPC\tgenotype-phenotype-assertions-IMPC.csv.gz\t2026-03-16\t2026-05-05",
    })

    assert rows == []


def test_phenotype_converter_skips_links_that_fail_legacy_2_of_3_support_rule():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {
        "IFX1": 101,
        "IFX2": 202,
    }

    rows = converter.phenotype_converter({
        "start_id": "NCBIGene:13489",
        "end_id": "MP:0005635",
        "protein_links": [
            {"protein_node_id": "IFX1", "support_sources": ["OMA"]},
            {"protein_node_id": "IFX2", "support_sources": ["HGNC", "Ensembl"]},
        ],
        "end_node": {"name": "decreased circulating bilirubin level"},
        "details": [{
            "source": "IMPC",
            "procedure_name": "Clinical Chemistry",
            "parameter_name": "Total bilirubin",
        }],
        "provenance": "IMPC\tgenotype-phenotype-assertions-IMPC.csv.gz\t2026-03-16\t2026-05-05",
    })

    assert rows == []


def test_phenotype_converter_maps_mgi_to_legacy_jax_ptype_without_nhprotein():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {
        "IFX1": 101,
        "IFX2": 202,
    }

    rows = converter.phenotype_converter({
        "start_id": "NCBIGene:13489",
        "end_id": "MP:0005385",
        "protein_links": [
            {"protein_node_id": "IFX1", "support_sources": ["OMA"]},
            {"protein_node_id": "IFX2", "support_sources": ["HGNC", "Ensembl"]},
        ],
        "end_node": {"name": "cardiovascular system phenotype"},
        "details": [{
            "source": "MGI",
            "source_id": "MGI:2166359",
        }],
        "provenance": "MGI\tHMD_HumanPhenotype.rpt\t2026-05-04\t2026-05-06",
    })

    assert len(rows) == 2
    assert {row.protein_id for row in rows} == {101, 202}
    assert all(row.ptype == "JAX/MGI Human Ortholog Phenotype" for row in rows)
    assert all(row.nhprotein_id is None for row in rows)
    assert all(row.term_id == "MP:0005385" for row in rows)
    assert all(row.term_name == "cardiovascular system phenotype" for row in rows)
    assert all(row.p_value is None for row in rows)
    assert all(row.procedure_name is None for row in rows)


def test_protein_mouse_phenotype_converter_maps_mgi_to_legacy_jax_ptype():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX1": 101}

    rows = converter.protein_mouse_phenotype_converter({
        "start_id": "IFX1",
        "end_id": "MP:0005385",
        "end_node": {"name": "cardiovascular system phenotype"},
        "details": [{
            "source": "MGI",
            "source_id": "MGI:2166359",
        }],
        "provenance": "MGI\tHMD_HumanPhenotype.rpt\t2026-05-04\t2026-05-06",
    })

    assert len(rows) == 1
    assert rows[0].ptype == "JAX/MGI Human Ortholog Phenotype"
    assert rows[0].protein_id == 101
    assert rows[0].nhprotein_id is None
    assert rows[0].term_id == "MP:0005385"
    assert rows[0].term_name == "cardiovascular system phenotype"


def test_ortholog_converter_emits_legacy_row_with_sorted_sources():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX1": 101}

    row = converter.ortholog_converter({
        "start_id": "IFX1",
        "species": "10090",
        "support_sources": ["EggNOG", "OMA", "Ensembl"],
        "source_db_ids": ["MGI:94924"],
        "ortholog_symbols": ["Drd2"],
        "ortholog_names": ["dopamine receptor D2"],
        "end_node": {
            "entrez_gene_id": "13489",
            "symbol": "Drd2",
            "name": "dopamine receptor D2",
        },
        "provenance": "HCOP\tlatest\t2025-05-16\t2026-05-05",
    })

    assert row.protein_id == 101
    assert row.taxid == 10090
    assert row.species == "Mus musculus"
    assert row.db_id == "MGI:94924"
    assert row.geneid == 13489
    assert row.symbol == "Drd2"
    assert row.name == "dopamine receptor D2"
    assert row.sources == "OMA, EggNOG"
    assert row.mod_url == "http://www.informatics.jax.org/marker/MGI:94924"


def test_ortholog_converter_skips_rows_failing_legacy_2_of_3_rule():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX1": 101}

    row = converter.ortholog_converter({
        "start_id": "IFX1",
        "species": "10090",
        "support_sources": ["OMA", "Ensembl"],
        "source_db_ids": ["MGI:94924"],
        "ortholog_symbols": ["Drd2"],
        "ortholog_names": ["dopamine receptor D2"],
        "end_node": {
            "entrez_gene_id": "13489",
            "symbol": "Drd2",
            "name": "dopamine receptor D2",
        },
        "provenance": "HCOP\tlatest\t2025-05-16\t2026-05-05",
    })

    assert row is None


def test_protein_converter_seeds_drgc_target_mapping_for_same_run():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    converter.protein_converter({
        "id": "IFX123",
        "xref": [],
        "gene_name": "SC6A5_HUMAN",
        "name": "Sodium- and chloride-dependent glycine transporter 2",
        "uniprot_id": "Q9Y345",
        "symbol": "SLC6A5",
        "sequence": "MSEQ",
        "dtoid": None,
        "dtoclass": None,
        "preferred_symbol": "SLC6A5",
        "provenance": "test",
    })

    row = converter.drgc_resource_converter({
        "id": "https://repo.metadatacenter.org/template-instances/example",
        "uniprot_id": "Q9Y345",
        "rssid": "https://repo.metadatacenter.org/template-instances/example",
        "resource_type": "Genetic Construct",
        "json": "{\"resource\": []}",
        "legacy_target_id": 999,
        "provenance": "Pharos 3.19\t3.19\t2024-02-15\tNone",
    })

    assert row.target_id == 123


def test_ligand_converter_assigns_integer_id_from_mapping():
    converter = TCRDOutputConverter()
    converter.id_mapping["ligand"] = {"LIGAND:EXISTING": 77}

    row = converter.ligand_converter({
        "id": "LIGAND:EXISTING",
        "name": "Existing ligand",
        "xref": [],
        "provenance": "test",
    })

    assert row.id == 77
    assert row.identifier == "LIGAND:EXISTING"


def test_entity_xrefs_preserve_full_ids_and_match_types_without_duplicates():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFXProtein:ONE": 11}
    rows = converter.protein_entity_xref_converter({
        "id": "IFXProtein:ONE",
        "xref": ["UniProtKB:P12345", "Symbol:OLD", "UniProtKB:P12345"],
        "xref_match_types": {
            "exact": ["UniProtKB:P12345"],
            "previous_symbol": ["Symbol:OLD"],
        },
    })
    assert [(row.entity_type, row.entity_id, row.source_id, row.match_type) for row in rows] == [
        ("protein", 11, "Symbol:OLD", "previous_symbol"),
        ("protein", 11, "UniProtKB:P12345", "exact"),
    ]


def test_entity_xref_dispatch_replaces_identity_xref_but_keeps_keyword_annotations():
    converter = TCRDOutputConverter()
    for model, method in [
        (Protein, 'protein_entity_xref_converter'),
        (Ligand, 'ligand_entity_xref_converter'),
        (Disease, 'disease_entity_xref_converter'),
    ]:
        converters = converter.get_object_converters(model)
        assert method in [convert.__name__ for convert in converters]
        assert not any(convert.__name__ == 'protein_xref_converter' for convert in converters)
    assert converter.protein_entity_xref_converter({'id': 'P', 'xref': ['UniProtKB:P1']})[0].__class__ is EntityXref
    assert converter.get_object_converters(ProteinKeywordEdge)[0].__name__ == 'keyword_xref_converter'
    assert converter.keyword_xref_converter({
        'start_id': 'P', 'end_node': {'source_id': 'KW-1'}, 'provenance': 'test'
    }).__class__ is Xref


def test_disease_entity_xrefs_only_include_exact_aliases():
    converter = TCRDOutputConverter()
    converter.id_mapping["ncats_disease"] = {"IFXDisease:ONE": 12}
    rows = converter.disease_entity_xref_converter({
        "id": "IFXDisease:ONE",
        "xref": ["IFXDisease:ONE", "MESH:D001943"],
        "xref_match_types": {"exact": ["IFXDisease:ONE", "MESH:D001943"]},
        "disease_harmonizer_payload": {"non_exact_xrefs": "MESH:D000001=broad"},
    })
    assert {row.source_id for row in rows} == {"IFXDisease:ONE", "MESH:D001943"}
    assert {row.match_type for row in rows} == {"exact"}


def test_ligand_entity_xrefs_include_ids_not_in_scalar_columns():
    converter = TCRDOutputConverter()
    converter.id_mapping["ligand"] = {"CHEBI:1": 13}
    rows = converter.ligand_entity_xref_converter({
        "id": "CHEBI:1",
        "xref": ["CHEBI:1", "CAS:123-45-6", "DRUGBANK:DB00001"],
    })
    assert {(row.entity_type, row.entity_id, row.source_id) for row in rows} == {
        ("ligand", 13, "CHEBI:1"),
        ("ligand", 13, "CAS:123-45-6"),
        ("ligand", 13, "DRUGBANK:DB00001"),
    }


def test_protein_converter_uses_short_graph_name_for_legacy_description():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    row = converter.protein_converter({
        "id": "IFX123",
        "xref": ["NCBIGene:1234", "ENSEMBL:ENSP00000123456"],
        "gene_name": "SC6A5_HUMAN",
        "name": "Sodium- and chloride-dependent glycine transporter 2",
        "description": "Long UniProt free text that should not land in legacy protein.description",
        "uniprot_id": "Q9Y345",
        "symbol": "SLC6A5",
        "sequence": "MSEQ",
        "dtoid": None,
        "dtoclass": None,
        "preferred_symbol": "SLC6A5",
        "provenance": "test",
    })

    assert row.name == "SC6A5_HUMAN"
    assert row.description == "Sodium- and chloride-dependent glycine transporter 2"


def test_tdl_info_converter_exports_uniprot_function_as_legacy_uniprot_function():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.tdl_info_converter({
        "id": "IFX123",
        "uniprot_function": "Long UniProt-derived descriptive text",
        "description": "Graph description fallback",
        "provenance": "test",
    })

    by_type = {row.itype: row for row in rows}
    assert by_type["UniProt Function"].protein_id == 123
    assert by_type["UniProt Function"].string_value == "Long UniProt-derived descriptive text"


def test_tdl_info_converter_exports_ncbi_gene_summary_as_legacy_summary():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.tdl_info_converter({
        "id": "IFX123",
        "ncbi_gene_summary": "NCBI-provided gene summary text",
        "provenance": "test",
    })

    by_type = {row.itype: row for row in rows}
    assert by_type["NCBI Gene Summary"].protein_id == 123
    assert by_type["NCBI Gene Summary"].string_value == "NCBI-provided gene summary text"


def test_ligand_edge_converter_uses_same_integer_ligand_id():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}
    converter.id_mapping["ligand"] = {"LIGAND:EXISTING": 77}

    rows = converter.ligand_edge_converter({
        "start_id": "IFX123",
        "end_id": "LIGAND:EXISTING",
        "details": [{
            "act_value": 1.5,
            "act_type": "IC50",
            "action_type": "inhibitor",
            "reference": "ref",
            "act_pmids": ["12345"],
        }],
        "provenance": "test",
    })

    assert len(rows) == 1
    assert rows[0].ncats_ligand_id == 77
    assert rows[0].target_id == 123


def test_pathway_converter_keeps_pwtype_without_lookup_table():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    row = converter.pathway_converter({
        "start_id": "IFX123",
        "end_node": {
            "type": "Reactome",
            "source_id": "R-HSA-199420",
            "name": "Generic Transcription Pathway",
            "url": "https://reactome.org/content/detail/R-HSA-199420",
        },
        "provenance": "Reactome\t95\t2025-11-27\t2026-03-23",
    })

    assert row.target_id == 123
    assert row.protein_id == 123
    assert row.pwtype == "Reactome"
    assert row.id_in_source == "R-HSA-199420"


def test_ppi_converter_emits_reciprocal_rows_with_max_score_and_legacy_stringdb_label():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {
        "IFX123": 123,
        "IFX456": 456,
    }

    rows = converter.ppi_converter({
        "start_id": "IFX456",
        "end_id": "IFX123",
        "sources": ["STRING\t12.0\t2023-05-16\t2026-04-16"],
        "score": [475, 477],
        "provenance": "STRING\t12.0\t2023-05-16\t2026-04-16",
    })

    assert len(rows) == 2
    assert {(row.protein_id, row.other_id) for row in rows} == {(123, 456), (456, 123)}
    assert all(row.ppitypes == "StringDB" for row in rows)
    assert all(row.score == 477 for row in rows)


def test_ppi_converter_joins_multiple_source_labels():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {
        "IFX123": 123,
        "IFX456": 456,
    }

    rows = converter.ppi_converter({
        "start_id": "IFX123",
        "end_id": "IFX456",
        "sources": [
            "Reactome\t95\t2025-11-27\t2026-03-23",
            "STRING\t12.0\t2023-05-16\t2026-04-16",
        ],
        "score": [800],
        "provenance": "STRING\t12.0\t2023-05-16\t2026-04-16",
    })

    assert len(rows) == 2
    assert all(row.ppitypes == "Reactome,StringDB" for row in rows)
    assert all(row.score == 800 for row in rows)


def test_ppi_converter_maps_reactome_pmids_and_interaction_type_downstream():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {
        "IFX123": 123,
        "IFX456": 456,
    }

    rows = converter.ppi_converter({
        "start_id": "IFX123",
        "end_id": "IFX456",
        "sources": ["Reactome\t96\t2026-03-24\t2026-04-24"],
        "pmids": [24243840, 11163199],
        "contexts": ["reactome:R-HSA-2428940"],
        "interaction_type": ["physical association"],
        "provenance": "Reactome\t96\t2026-03-24\t2026-04-24",
    })

    assert len(rows) == 2
    assert all(row.evidence == "24243840|11163199" for row in rows)
    assert all(row.interaction_type == "physical association" for row in rows)


def test_virus_converter_maps_graph_node_to_legacy_virus_row():
    converter = TCRDOutputConverter()

    row = converter.virus_converter({
        "id": "NCBITaxon:1006008",
        "nucleic1": "RNA",
        "nucleic2": "ssRNA(+)",
        "order": "Picornavirales",
        "family": "Picornaviridae",
        "subfamily": None,
        "genus": "Enterovirus",
        "species": "Enterovirus A",
        "name": "Human enterovirus 71 HZ08/Hangzhou/2008",
    })

    assert row.virusTaxid == "1006008"
    assert row.family == "Picornaviridae"
    assert row.name == "Human enterovirus 71 HZ08/Hangzhou/2008"


def test_viral_protein_converter_uses_preloaded_graph_edge_map_for_virus_id():
    converter = TCRDOutputConverter()
    converter.id_mapping["viral_protein_to_virus"] = {
        "PHIPSTER.ViralProtein:15": "1006008",
    }

    row = converter.viral_protein_converter({
        "id": "PHIPSTER.ViralProtein:15",
        "name": "full_polyprotein 1..2193",
        "ncbi": "AEB71507.1",
    })

    assert row.id == 15
    assert row.name == "full_polyprotein 1..2193"
    assert row.ncbi == "AEB71507.1"
    assert row.virus_id == "1006008"


def test_viral_ppi_converter_round_trips_list_pdb_ids_to_legacy_scalar_column():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"UniProtKB:P19320": 123}

    rows = converter.viral_ppi_converter({
        "start_id": "UniProtKB:P19320",
        "end_id": "PHIPSTER.ViralProtein:15",
        "details": [
            {
                "source": "P-HIPSTer",
                "source_protein_id": "UniProtKB:P19320",
                "final_lr": 128.345190881009,
                "pdb_ids": ["1z7z:3I(0.09,0.19,19/20,12,14)"],
                "high_confidence": False,
            },
            {
                "source": "P-HIPSTer",
                "source_protein_id": "UniProtKB:P19320",
                "final_lr": 408.948397259193,
                "pdb_ids": [],
                "high_confidence": True,
            },
        ],
        "provenance": "P-HIPSTer\toriginal paper\t2019-09-05\t2026-05-13",
    })

    assert len(rows) == 2
    assert rows[0].viral_protein_id == 15
    assert rows[0].protein_id == 123
    assert rows[0].dataSource == "P-HIPSTer"
    assert rows[0].finalLR == 128.345190881009
    assert rows[0].pdbIDs == "1z7z:3I(0.09,0.19,19/20,12,14)"
    assert rows[0].highConfidence is False
    assert rows[1].finalLR == 408.948397259193
    assert rows[1].pdbIDs is None
    assert rows[1].highConfidence is True


def test_gtex_converter_branches_gtex_details_from_shared_expression_edge():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    rows = converter.gtex_converter({
        "start_id": "IFX123",
        "details": [
            {
                "source": "GTEx",
                "tissue": "Liver",
                "source_tissue_id": "UBERON:0002107",
                "number_value": 10.0,
                "source_rank": 0.7,
            },
            {
                "source": "GTEx",
                "tissue": "Liver",
                "source_tissue_id": "UBERON:0002107",
                "sex": "male",
                "number_value": 11.0,
                "source_rank": 0.8,
            },
            {
                "source": "GTEx",
                "tissue": "Liver",
                "source_tissue_id": "UBERON:0002107",
                "sex": "female",
                "number_value": 9.0,
                "source_rank": 0.6,
            },
            {
                "source": "HPM Protein",
                "tissue": "Liver",
                "source_tissue_id": "UBERON:0002107",
                "number_value": 1.0,
            },
        ],
        "provenance": "GTEx\tv1\t2026-01-01\t2026-01-02",
    })

    assert len(rows) == 1
    row = rows[0]
    assert row.protein_id == 123
    assert row.tissue == "Liver"
    assert row.uberon_id == "UBERON:0002107"
    assert row.tpm == 10.0
    assert row.tpm_male == 11.0
    assert row.tpm_female == 9.0


def test_expression_converter_uses_resolved_tissue_id_for_uberon_id():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}
    converter.id_mapping["tissue"] = {"Brain": 85}

    rows = converter.expression_converter({
        "start_id": "IFX123",
        "end_id": "UBERON:0000955",
        "details": [
            {
                "source": "JensenLab TISSUES",
                "tissue": "Brain",
                "source_id": "ENSP00000000001",
                "source_tissue_id": "BTO:0000142",
                "number_value": 3.5,
                "expressed": True,
            },
        ],
        "provenance": "JensenLab TISSUES\tNone\t2026-01-01\t2026-01-02",
    })

    assert len(rows) == 1
    assert rows[0].etype == "JensenLab TISSUES"
    assert rows[0].tissue_id == 85
    assert rows[0].uberon_id == "UBERON:0000955"


def test_mondo_xref_converter_marks_exact_matches():
    converter = TCRDOutputConverter()

    rows = converter.mondo_xref_converter({
        "id": "MONDO:0000001",
        "mondo_xrefs": ["DOID:4", "ICD9:799.9"],
        "exact_matches": ["DOID:4"],
        "provenance": "Mondo\tv2026-03-03\t2026-03-03\t2026-03-10",
    })

    assert len(rows) == 2
    by_xref = {row.xref: row for row in rows}
    assert by_xref["DOID:4"].db == "DOID"
    assert by_xref["DOID:4"].value == "4"
    assert by_xref["DOID:4"].equiv_to is True
    assert by_xref["ICD9:799.9"].equiv_to is False


def test_keyword_xref_converter_writes_keyword_id_to_value_and_label_to_xtra():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    row = converter.keyword_xref_converter({
        "start_id": "IFX123",
        "end_node": {
            "id": "keyword:uniprot:technical term:reference proteome",
            "source_id": "KW-1185",
            "value": "Reference proteome",
            "category": "Technical term",
        },
        "provenance": "UniProt\t2026_01\t2026-01-28\t2026-03-16",
    })

    assert row.protein_id == 123
    assert row.xtype == "UniProt Keyword"
    assert row.value == "KW-1185"
    assert row.xtra == "Reference proteome"


def test_do_parent_uses_composite_primary_key():
    pk_columns = tuple(column.name for column in DOParent.__table__.primary_key.columns)

    assert pk_columns == ("doid", "parent_id")


def test_p2p_converter_dedupes_across_calls():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}

    obj = {
        "start_id": "IFX123",
        "gene_id": 6857,
        "pmids": ["20222955", "20222955"],
        "provenance": "Pharos 4.0 CSV\t1.0\tNone\t2024-10-03",
    }

    first_rows = converter.p2p_converter(obj)
    second_rows = converter.p2p_converter(obj)

    assert len(first_rows) == 1
    assert len(second_rows) == 0
    assert first_rows[0].protein_id == 123
    assert str(first_rows[0].pubmed_id) == "20222955"
    assert first_rows[0].gene_id == 6857
    assert first_rows[0].source == "NCBI"


def test_ncats_disease_converter_only_keeps_loaded_mondo_ids():
    converter = TCRDOutputConverter()

    missing = converter.ncats_disease_converter({
        "id": "MONDO:0957221",
        "name": "Spastic paraplegia 70, autosomal recessive",
        "provenance": "Mondo\tv2026-03-03\t2026-03-03\t2026-03-10",
    })
    assert missing.mondoid is None

    converter._known_mondo_ids.add("MONDO:0000001")
    known = converter.ncats_disease_converter({
        "id": "MONDO:0000001",
        "name": "disease",
        "provenance": "Mondo\tv2026-03-03\t2026-03-03\t2026-03-10",
    })
    assert known.mondoid == "MONDO:0000001"


def test_ncats_disease_converter_maps_rare_disease_to_gard_rare():
    converter = TCRDOutputConverter()

    row = converter.ncats_disease_converter({
        "id": "MONDO:0000001",
        "name": "rare disease",
        "rare_disease": True,
        "provenance": "RDAS\trdas-diseases-graphql\tNone\t2026-05-27",
    })

    assert row.gard_rare is True


def test_disease_converter_prefers_detail_source_id_and_keeps_resolved_mondoid():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}
    converter.id_mapping["ncats_disease"] = {"MONDO:0000001": 456}
    converter._known_mondo_ids.add("MONDO:0000001")

    rows = converter.disease_converter({
        "start_id": "IFX123",
        "end_id": "MONDO:0000001",
        "end_node": {"name": "Example disease"},
        "details": [{
            "source": "JensenLab Text Mining",
            "source_id": "DOID:1234",
            "confidence": 3.2,
            "zscore": 6.5,
            "url": "https://example.org/jensen",
        }],
        "provenance": "JensenLab DISEASES\tNone\t2026-03-18\t2026-04-02",
    })

    assert len(rows) == 1
    row = rows[0]
    assert row.protein_id == 123
    assert row.dtype == "JensenLab Text Mining"
    assert row.did == "DOID:1234"
    assert row.mondoid == "MONDO:0000001"
    assert row.conf == 3.2
    assert row.zscore == 6.5
    assert row.reference == "https://example.org/jensen"


def test_disease_converter_carries_original_sources_into_mysql_source():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"UniProtKB:P11111": 123}

    rows = converter.disease_converter({
        "start_id": "UniProtKB:P11111",
        "end_id": "DOID:1",
        "end_node": {"name": "Disease One"},
        "details": [{
            "source": "eRAM",
            "source_id": "DOID:1",
            "original_sources": ["CTD_human", "ORPHANET", "CLINVAR"],
        }],
        "provenance": "Pharos 3.19\t3.19\t2024-02-15\tNone",
    })

    assert len(rows) == 1
    row = rows[0]
    assert row.dtype == "eRAM"
    assert row.did == "DOID:1"
    assert row.evidence is None
    assert row.source == "CTD_human|ORPHANET|CLINVAR"


def test_ncats_d2da_converter_keys_links_off_resolved_disease_id():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFX123": 123}
    converter.id_mapping["ncats_disease"] = {"MONDO:0000001": 456}
    converter._known_mondo_ids.add("MONDO:0000001")

    obj = {
        "start_id": "IFX123",
        "end_id": "MONDO:0000001",
        "end_node": {"name": "Example disease"},
        "details": [{
            "source": "JensenLab Text Mining",
            "source_id": "DOID:1234",
            "confidence": 3.2,
            "zscore": 6.5,
            "url": "https://example.org/jensen",
        }],
        "provenance": "JensenLab DISEASES\tNone\t2026-03-18\t2026-04-02",
    }

    disease_rows = converter.disease_converter(obj)
    link_rows = converter.ncats_d2da_converter(obj)

    assert len(disease_rows) == 1
    assert len(link_rows) == 1
    assert link_rows[0].ncats_disease_id == 456
    assert link_rows[0].disease_assoc_id == disease_rows[0].id
    assert link_rows[0].direct == 1


def test_promoted_gene_disease_association_exports_to_named_disease():
    detail = DiseaseAssociationDetail(
        source="ClinGen",
        source_id="MONDO:0000001",
        source_subject_id="HGNC:1234",
        subject_type="gene",
    )
    gene_edge = GeneDiseaseEdge(
        start_node=Gene(id="IFXGene:SOURCE"),
        end_node=Disease(id="IFXDisease:DEST", name="Example disease"),
        details=[detail],
        gene_assessments=[{"gene_id": "IFXGene:SOURCE", "tier": "curated"}],
    )
    promoted = InputAdapter._canonicalize_relationship_class(
        gene_edge,
        Protein(id="IFXProtein:TARGET"),
        gene_edge.end_node,
    )

    assert isinstance(promoted, ProteinDiseaseEdge)
    assert promoted.details == [detail]
    assert promoted.gene_assessments == gene_edge.gene_assessments

    converter = TCRDOutputConverter()
    converter.id_mapping = {
        "protein": {"IFXProtein:TARGET": 123},
        "ncats_disease": {"IFXDisease:DEST": 456},
    }
    converter._disease_mondoid_by_id["IFXDisease:DEST"] = "MONDO:0000001"
    row = {
        "start_id": promoted.start_node.id,
        "end_id": promoted.end_node.id,
        "end_node": {"name": promoted.end_node.name},
        "details": [
            detail.to_dict(),
            DiseaseAssociationDetail(
                source="ClinGen",
                source_id="MONDO:0000001",
                source_subject_id="HGNC:1234",
                subject_type="gene",
                negated=True,
            ).to_dict(),
        ],
        "gene_assessments": promoted.gene_assessments,
        "provenance": "DISEASES\t2.3.5",
    }
    associations = converter.disease_converter(row)
    links = converter.ncats_d2da_converter(row)
    assert len(associations) == len(links) == 1
    association = associations[0]
    link = links[0]

    assert association.protein_id == 123
    assert association.dtype == "ClinGen"
    assert association.did == "MONDO:0000001"
    assert link.ncats_disease_id == 456
    assert link.disease_assoc_id == association.id
    assert link.direct == 1

    row["details"] = [row["details"][1]]
    assert converter.disease_type_converter(row) is None
    assert converter.disease_converter(row) == []
    assert converter.ncats_d2da_converter(row) == []


def test_affiliate_converter_uses_provider_name_for_display_name():
    converter = TCRDOutputConverter()

    row = converter.affiliate_converter({
        "source": "TIGA",
        "name": "Target Illumination GWAS Analytics (TIGA)",
        "description": "TIGA scores and ranks GWAS discovered associations.",
        "link_count": 18721,
    })

    assert row.source == "TIGA"
    assert row.display_name == "Target Illumination GWAS Analytics (TIGA)"
    assert row.description == "TIGA scores and ranks GWAS discovered associations."
    assert row.link_count == 18721


def test_extlink_converter_picks_lexicographically_first_source_id_detail():
    converter = TCRDOutputConverter()
    converter.id_mapping["protein"] = {"IFXProtein:1": 123}

    row = converter.extlink_converter({
        "start_id": "IFXProtein:1",
        "source": "TIGA",
        "details": [
            {
                "source_id": "ENSG000002",
                "url": "https://unmtid-shinyapps.net/shiny/tiga/?gene=ENSG000002",
            },
            {
                "source_id": "ENSG000001",
                "url": "https://unmtid-shinyapps.net/shiny/tiga/?gene=ENSG000001",
            },
        ],
    })

    assert row.protein_id == 123
    assert row.source == "TIGA"
    assert row.url == "https://unmtid-shinyapps.net/shiny/tiga/?gene=ENSG000001"
