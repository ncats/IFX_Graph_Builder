-- RaMP 3.0.7 application schema; redundant reaction_protein2met omitted.
CREATE TABLE "analyte" (
	"rampId" VARCHAR(30) NOT NULL  ,
	"type" VARCHAR(30) NULL  , common_name TEXT,
	PRIMARY KEY ("rampId")
);

CREATE TABLE "analytehasontology"
(
    rampCompoundId VARCHAR(30)
        constraint analytehasontology_analyte_rampId_fk
            references analyte,
    rampOntologyId VARCHAR(30)
        constraint analytehasontology_ontology_rampOntologyId_fk
            references ontology,
    constraint analytehasontology_pk
        primary key (rampCompoundId, rampOntologyId)
);

CREATE TABLE "analytehaspathway"
(
    rampId        VARCHAR(30)
        constraint analytehaspathway_analyte_rampId_fk
            references analyte,
    pathwayRampId VARCHAR(30)
        constraint analytehaspathway_pathway_pathwayRampId_fk
            references pathway,
    pathwaySource VARCHAR(30)
);

CREATE TABLE "analytesynonym"
(
    Synonym        varchar(500) collate NOCASE,
    rampId         varchar(30)
        constraint analytesynonym_analyte_rampId_fk
            references analyte,
    geneOrCompound varchar(30),
    source         varchar(30)
);

CREATE TABLE "catalyzed"
(
    rampCompoundId VARCHAR(30)
        constraint catalyzed_analyte_rampId_fk
            references analyte,
    rampGeneId     VARCHAR(30)
        constraint catalyzed_analyte_rampId_fk2
            references analyte,
    proteinType    VARCHAR(32),
    constraint catalyzed_pk
        primary key (rampCompoundId, rampGeneId)
);

CREATE TABLE "chem_props"
(
    ramp_id          VARCHAR(30) not null
        constraint chem_props_analyte_rampId_fk
            references analyte,
    chem_data_source VARCHAR(32),
    chem_source_id   VARCHAR(45),
    iso_smiles       VARCHAR(256),
    inchi_key_prefix VARCHAR(32),
    inchi_key        VARCHAR(32),
    inchi            VARCHAR(4096),
    mw               FLOAT,
    monoisotop_mass  FLOAT,
    common_name      VARCHAR(1024),
    mol_formula      VARCHAR(64)
);

CREATE TABLE "db_version" (
	"ramp_version" VARCHAR(20) NOT NULL  ,
	"load_timestamp" DATETIME NOT NULL DEFAULT 'CURRENT_TIMESTAMP' ,
	"version_notes" VARCHAR(256) NULL  ,
	"met_intersects_json" VARCHAR(10000) NULL  ,
	"gene_intersects_json" VARCHAR(10000) NULL  ,
	"met_intersects_json_pw_mapped" VARCHAR(10000) NULL  ,
	"gene_intersects_json_pw_mapped" VARCHAR(10000) NULL  ,
	"db_sql_url" VARCHAR(256) NULL
);

CREATE TABLE "entity_status_info" (
	"status_category" VARCHAR(64) NOT NULL  ,
	"entity_source_id" VARCHAR(32) NOT NULL  ,
	"entity_source_name" VARCHAR(45) NOT NULL  ,
	"entity_count" INTEGER NOT NULL
);

CREATE TABLE "metabolite_class"
(
    ramp_id          VARCHAR(32)  not null
        constraint metabolite_class_analyte_rampId_fk
            references analyte,
    class_source_id  VARCHAR(32)  not null,
    class_level_name VARCHAR(128) not null,
    class_name       VARCHAR(128) not null,
    source           VARCHAR(32)  not null,
    constraint metabolite_class_pk
        primary key (class_source_id, ramp_id, class_level_name)
);

CREATE TABLE "ontology"
(
    rampOntologyId   varchar(30)
        constraint ontology_pk
            primary key,
    commonName       varchar(64) collate NOCASE,
    HMDBOntologyType varchar(30),
    metCount         int
);

CREATE TABLE pathway (
 pathwayRampId varchar(30) primary key,
 sourceId varchar(30),
 type varchar(30),
 pathwayCategory varchar(30),
 pathwayName varchar(250) COLLATE NOCASE);

CREATE TABLE pathway_duplicates
(
    pathwayRampId1 varchar(30) not null,
    pathwayRampId2 varchar(30) not null
);

CREATE TABLE "pathway_similarity"
(
    pathwayRampId    varchar(30) not null
        constraint pathway_similarity_pk
            primary key,
    analyte_blob     blob,
    metabolite_blob  blob,
    gene_blob        blob,
    metabolite_count integer default 0,
    gene_count       integer default 0
);

CREATE TABLE "reaction" (
	"ramp_rxn_id" VARCHAR(16) NOT NULL  ,
	"rxn_source_id" VARCHAR(16) NOT NULL  ,
	"status" INTEGER NOT NULL  ,
	"is_transport" INTEGER NOT NULL  ,
	"direction" VARCHAR(8) NOT NULL  ,
	"label" VARCHAR(256) NOT NULL  ,
	"equation" VARCHAR(256) NOT NULL  ,
	"html_equation" VARCHAR(256) NOT NULL  ,
	"ec_num" VARCHAR(256) NULL  ,
	"has_human_prot" INTEGER NOT NULL  ,
	"only_human_mets" INTEGER NOT NULL  ,
	PRIMARY KEY ("ramp_rxn_id")
);

CREATE TABLE "reaction2met"
(
    ramp_rxn_id       VARCHAR(16)         not null
        constraint reaction2met_reaction_ramp_rxn_id_fk
            references reaction,
    rxn_source_id     VARCHAR(16)         not null,
    ramp_cmpd_id      VARCHAR(16)         not null
        constraint reaction2met_analyte_rampId_fk
            references analyte,
    substrate_product INTEGER             not null,
    met_source_id     VARCHAR(32)         not null,
    met_name          VARCHAR(256),
    is_cofactor       INTEGER default '0' not null
);

CREATE TABLE "reaction2protein"
(
    ramp_rxn_id   VARCHAR(16)       not null
        constraint reaction2protein_reaction_ramp_rxn_id_fk
            references reaction,
    rxn_source_id VARCHAR(16)       not null,
    ramp_gene_id  VARCHAR(16)       not null
        constraint reaction2protein_analyte_rampId_fk
            references analyte,
    uniprot       VARCHAR(16)       not null,
    protein_name  VARCHAR(16)       not null,
    is_reviewed   integer default 0 not null
);

CREATE TABLE "reaction_ec_class"
(
    ramp_rxn_id         varchar(16)  not null
        constraint reaction_ec_class_reaction_ramp_rxn_id_fk
            references reaction,
    rxn_source_id       varchar(16)  not null,
    rxn_class_ec        varchar(16)  not null,
    ec_level            int          not null,
    rxn_class           varchar(256) not null,
    rxn_class_hierarchy varchar(512) not null
);

CREATE TABLE "source"
(
    sourceId           VARCHAR(30)         not null,
    rampId             VARCHAR(30)
        constraint source_analyte_rampId_fk
            references analyte,
    IDtype             VARCHAR(30),
    geneOrCompound     VARCHAR(30),
    commonName         VARCHAR(256),
    priorityHMDBStatus VARCHAR(32),
    dataSource         VARCHAR(32),
    pathwayCount       INTEGER default '0' not null
);

CREATE TABLE "version_info" (
	"ramp_db_version" VARCHAR(16) NOT NULL  ,
	"db_mod_date" DATE NOT NULL  ,
	"status" VARCHAR(16) NOT NULL  ,
	"data_source_id" VARCHAR(32) NOT NULL  ,
	"data_source_name" VARCHAR(128) NOT NULL  ,
	"data_source_url" VARCHAR(128) NOT NULL  ,
	"data_source_version" VARCHAR(128) NOT NULL
);

CREATE INDEX ahp_RampID_IDX
    on analytehaspathway (rampId);

CREATE INDEX ahp_path_source_IDX
    on analytehaspathway (pathwaySource);

CREATE INDEX analyte_ont_id_idx
    on analytehasontology (rampOntologyId);

CREATE INDEX analyte_ont_ramp_id_idx
    on analytehasontology (rampCompoundId);

CREATE INDEX "analyte_rampId_RampID_IDX" ON "analyte" ("rampId");

CREATE INDEX analytesynonym_Synonym_index
    on analytesynonym (Synonym);

CREATE INDEX catal_comp_idx
    on catalyzed (rampCompoundId);

CREATE INDEX catal_gene_idx
    on catalyzed (rampGeneId);

CREATE INDEX class_name_metclass_idx
    on metabolite_class (class_name);

CREATE INDEX class_source_id_metclass_idx
    on metabolite_class (class_source_id);

CREATE INDEX class_source_metclass_idx
    on metabolite_class (source);

CREATE INDEX "data_source_index" ON "version_info" ("data_source_id");

CREATE INDEX ec_class_ramp_rxn_id
    on reaction_ec_class (ramp_rxn_id);

CREATE INDEX ec_class_rxn_source_idx
    on reaction_ec_class (rxn_source_id);

CREATE INDEX ec_level_idx
    on reaction_ec_class (ec_level);

CREATE INDEX idx_is_reviewed_rxn2prot
    on reaction2protein (is_reviewed);

CREATE INDEX idx_source_IDtype
    on source (IDtype);

CREATE INDEX idx_source_geneOrCompound
    on source (geneOrCompound);

CREATE INDEX inchi_key_idx
    on chem_props (inchi_key);

CREATE INDEX inchi_key_prefix_idx
    on chem_props (inchi_key_prefix);

CREATE INDEX pathwayRampID_IDX
    on analytehaspathway (pathwayRampId);

CREATE INDEX prop_source_idx
    on chem_props (chem_data_source);

CREATE INDEX ramp_id_idx
    on chem_props (ramp_id);

CREATE INDEX ramp_id_metclass_idx
    on metabolite_class (ramp_id);

CREATE INDEX "reaction_ec_num_idx" ON "reaction" ("ec_num");

CREATE INDEX "reaction_has_human_prot_idx" ON "reaction" ("has_human_prot");

CREATE INDEX "reaction_src_id_idx" ON "reaction" ("rxn_source_id");

CREATE INDEX rxn2met_iscofactor_idx
    on reaction2met (is_cofactor);

CREATE INDEX rxn2met_met_ramp_id_idx
    on reaction2met (ramp_cmpd_id);

CREATE INDEX rxn2met_met_source_id_idx
    on reaction2met (met_source_id);

CREATE INDEX rxn2met_rxn_ramp_id_idx
    on reaction2met (ramp_rxn_id);

CREATE INDEX rxn2met_rxn_source_id_idx
    on reaction2met (rxn_source_id);

CREATE INDEX rxn2met_subs_prod_idx
    on reaction2met (substrate_product);

CREATE INDEX rxn2prot_ramp_gene_id_idx
    on reaction2protein (ramp_gene_id);

CREATE INDEX rxn2prot_source_id_idx
    on reaction2protein (rxn_source_id);

CREATE INDEX rxn2prot_uniprot_idx
    on reaction2protein (uniprot);

CREATE INDEX rxn_prot_ramp_id_idx
    on reaction2protein (ramp_rxn_id);

CREATE INDEX source_RampID_IDX
    on source (rampId);

CREATE INDEX source_comName_RampID_IDX
    on source (commonName);

CREATE INDEX source_datasrc_IDX
    on source (dataSource);

CREATE INDEX source_pathCount_IDX
    on source (pathwayCount);

CREATE INDEX source_sid_RampID_IDX
    on source (sourceId);

CREATE INDEX "status_index" ON "version_info" ("status");

ALTER TABLE version_info ADD COLUMN data_source_snapshot_ids TEXT;
CREATE TABLE ramp_export_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
