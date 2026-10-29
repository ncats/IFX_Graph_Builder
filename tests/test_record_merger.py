from src.shared.record_merger import RecordMerger, FieldConflictBehavior


def test_gene_assessments_merge_on_shared_protein_edge():
    records = [
        {
            "start_id": "IFXProtein:P", "end_id": "IFXDisease:D",
            "gene_assessments": [{"gene_id": "IFXGene:A", "tier": "refuted", "score": "0", "origins": ["ClinGen"]}],
            "entity_resolution": "gene A", "provenance": "DISEASES 2.3.3",
        },
        {
            "start_id": "IFXProtein:P", "end_id": "IFXDisease:D",
            "gene_assessments": [{"gene_id": "IFXGene:B", "tier": "supported", "score": "4", "origins": ["OMIM"]}],
            "entity_resolution": "gene B", "provenance": "DISEASES 2.3.3",
        },
    ]
    merged = RecordMerger().merge_records(records, {}, nodes_or_edges="edges")
    assert len(merged) == 1
    assert [value["gene_id"] for value in merged[0]["gene_assessments"]] == ["IFXGene:A", "IFXGene:B"]

def get_existing_records():
    return {
        "1": {
            "id": "1",
            "name": "Aleece",
            "old_field": "value2",
            "old_list_field": ["old_value"],
            "resolved_ids": ["res2"],
            "xref": ["existing_xref"],
            "creation": "source2",
            "updates": ["update1"]
        }
    }

def get_merging_records():
    return [
        {
            "id": "1",
            "name": "Alice",
            "old_list_field": ["another_value"],
            "new_field": "value1",
            "new_list_field": ["new_value"],
            "empty_list_field": [],
            "entity_resolution": "res1",
            "xref": ["xref1"],
            "provenance": "source1"
        }
    ]

def test_basic_record_merging():

    for field_conflict_behavior in [FieldConflictBehavior.KeepFirst, FieldConflictBehavior.KeepLast]:
        # Arrange
        existing_record_map = get_existing_records()
        records = get_merging_records()
        merger = RecordMerger(field_conflict_behavior=field_conflict_behavior)

        # Act
        merged_records = merger.merge_records(records, existing_record_map)

        # Assert

        assert len(merged_records) == 1
        assert merged_records[0]["id"] == "1"
        assert merged_records[0]["name"] == "Alice" if field_conflict_behavior == FieldConflictBehavior.KeepLast else "Aleece"
        assert merged_records[0]["old_field"] == "value2"
        assert sorted(merged_records[0]["old_list_field"]) == sorted(["old_value", "another_value"])
        assert merged_records[0]["new_field"] == "value1"
        assert merged_records[0]["new_list_field"] == ["new_value"]
        assert "empty_list_field" not in merged_records[0]
        assert sorted(merged_records[0]["resolved_ids"]) == sorted(["res2", "res1"])
        assert merged_records[0]["xref"] == ["existing_xref"]
        assert merged_records[0]['creation'] == "source2"
        assert sorted(merged_records[0]["updates"]) == sorted([
            'update1',
            f'name\tAleece\tAlice\tsource1\t{field_conflict_behavior.value}',
            f'new_field\tNULL\tvalue1\tsource1\t{field_conflict_behavior.value}',
            'new_list_field\tNULL\tadding ["new_value"]\tsource1',
            'old_list_field\t1 entries already there\tadding ["another_value"]\tsource1'
        ])

def test_new_record_merging():

    # Arrange
    existing_record_map = {}
    records = get_merging_records()
    merger = RecordMerger(field_conflict_behavior=FieldConflictBehavior.KeepLast)

    # Act
    merged_records = merger.merge_records(records, existing_record_map)

    assert len(merged_records) == 1
    assert merged_records[0]["id"] == "1"
    assert merged_records[0]["name"] == "Alice"
    assert "old_field" not in merged_records[0]
    assert merged_records[0]["old_list_field"] == ["another_value"]
    assert merged_records[0]["new_field"] == "value1"
    assert merged_records[0]["new_list_field"] == ["new_value"]
    assert "empty_list_field" not in merged_records[0]
    assert merged_records[0]["resolved_ids"] == ["res1"]
    assert merged_records[0]["xref"] == ["xref1"]
    assert merged_records[0]["creation"] == "source1"
    assert "updates" not in merged_records[0]

def test_null_merge():
    # Arrange
    existing_record_map = get_existing_records()
    records = [
        {"id": "2", "name": "Alice", "first_has": "val", "first_hasnt": None, "first_has_list": [1], "first_hasnt_list": [], "provenance": "source1", "entity_resolution": "res1"},
        {"id": "3", "name": "Alice", "first_has": None, "first_hasnt": "val", "first_has_list": [], "first_hasnt_list": [2], "provenance": "source1", "entity_resolution": "res2"}
    ]
    merger = RecordMerger(field_conflict_behavior=FieldConflictBehavior.KeepLast)

    # Act
    merged_records = merger.merge_records(records, existing_record_map)

    # Assert
    assert len(merged_records) == 3
    assert merged_records[1]["id"] == "2"
    assert merged_records[1]["first_has"] == "val"
    assert "first_hasnt" not in merged_records[1]
    assert merged_records[1]["first_has_list"] == [1]
    assert "first_hasnt_list" not in merged_records[1]

    assert merged_records[2]["first_hasnt"] == "val"
    assert "first_has" not in merged_records[2]
    assert merged_records[2]["first_hasnt_list"] == [2]
    assert "first_has_list" not in merged_records[2]



def test_double_merge():
    # Arrange
    existing_record_map = get_existing_records()
    records = [
        {"id": "1", "name": "Alice", "new_field": "value1", "entity_resolution": "res1", "provenance": "source1", "xref": ["xref1"], "old_list_field": ["another_value"]},
        {"id": "1", "name": "Aleece", "new_field": "value1", "entity_resolution": "res2", "provenance": "source2", "xref": ["xref2"], "old_list_field": ["yet_another_value"]}
    ]
    merger = RecordMerger(field_conflict_behavior=FieldConflictBehavior.KeepLast)

    # Act
    merged_records = merger.merge_records(records, existing_record_map)

    # Assert
    assert len(merged_records) == 1
    assert merged_records[0]["id"] == "1"
    assert merged_records[0]["name"] == "Aleece"
    assert merged_records[0]["old_field"] == "value2"
    assert sorted(merged_records[0]["old_list_field"]) == sorted(["old_value", "another_value", "yet_another_value"])
    assert merged_records[0]["new_field"] == "value1"
    assert "empty_list_field" not in merged_records[0]
    assert merged_records[0]["xref"] == ["existing_xref"]
    assert merged_records[0]['creation'] == "source2"
    assert sorted(merged_records[0]["resolved_ids"]) == sorted(["res2", "res1"])
    assert merged_records[0]["updates"] == [
        'update1',
        'name\tAleece\tAlice\tsource1\tKeepLast',
        'new_field\tNULL\tvalue1\tsource1\tKeepLast',
        'old_list_field\t1 entries already there\tadding ["another_value"]\tsource1',
        'name\tAlice\tAleece\tsource2\tKeepLast',
        'old_list_field\t2 entries already there\tadding ["yet_another_value"]\tsource2',
    ]


def test_list_update_audit_caps_added_entries_preview_at_five():
    existing_record_map = {
        "1": {
            "id": "1",
            "sources": ["existing-source"],
            "resolved_ids": ["res-existing"],
            "creation": "source-existing",
            "updates": []
        }
    }
    records = [{
        "id": "1",
        "sources": [
            "existing-source",
            "source-1",
            "source-2",
            "source-3",
            "source-4",
            "source-5",
            "source-6",
        ],
        "entity_resolution": "res-new",
        "provenance": "source-new"
    }]
    merger = RecordMerger(field_conflict_behavior=FieldConflictBehavior.KeepLast)

    merged_records = merger.merge_records(records, existing_record_map)

    assert merged_records[0]["updates"] == [
        'sources\t1 entries already there\tadding ["source-1", "source-2", "source-3", "source-4", "source-5"] (+1 more)\tsource-new'
    ]


def test_record_merger_dedupes_list_of_dict_patent_family_mentions():
    existing_record_map = {
        "IFX1": {
            "id": "IFX1",
            "patent_family_mentions": [
                "2020:1001",
                "2021:1002",
            ],
            "resolved_ids": ["res-existing"],
            "creation": "source-existing",
            "updates": [],
        }
    }
    records = [{
        "id": "IFX1",
        "patent_family_mentions": [
            "2020:1001",
            "2021:1003",
        ],
        "entity_resolution": "res-new",
        "provenance": "source-new",
    }]
    merger = RecordMerger(field_conflict_behavior=FieldConflictBehavior.KeepLast)

    merged_records = merger.merge_records(records, existing_record_map)

    assert merged_records[0]["patent_family_mentions"] == [
        "2020:1001",
        "2021:1002",
        "2021:1003",
    ]
