from sqlalchemy import Column, String, Integer, Text, Enum, SmallInteger, ForeignKey, Date, DateTime, Float, Index, Boolean, \
    DECIMAL, UniqueConstraint, Double, text
from sqlalchemy.dialects.mysql import ENUM
from sqlalchemy.orm import declarative_base, declared_attr
from enum import Enum as PyEnum


class TDL(PyEnum):
    TclinPlus = 'Tclin+'
    Tclin = 'Tclin'
    TchemPlus = 'Tchem+'
    Tchem = 'Tchem'
    Tbio = 'Tbio'
    Tgray = 'Tgray'
    Tdark = 'Tdark'

class GoType(PyEnum):
    Component = 'Component'
    Function = 'Function'
    Process = 'Process'

class BaseMixin:
    @declared_attr
    def provenance(cls):
        return Column(Text)

Base = declarative_base(cls=BaseMixin)

class HgramCDF(Base):
    __tablename__ = "hgram_cdf"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    type = Column(String(255), ForeignKey("gene_attribute_type.name"), nullable=False)
    attr_count = Column(Integer, nullable=False)
    attr_cdf = Column(DECIMAL(17, 16), nullable=False)

    __table_args__ = (
        Index("hgram_cdf_idx1", "protein_id"),
        Index("hgram_cdf_idx2", "type"),
    )

class Protein(Base):
    __tablename__ = "protein"

    id = Column(Integer, primary_key=True, nullable=False, autoincrement=False)
    ifx_id = Column(String(18), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=False)
    uniprot = Column(String(20))
    up_version = Column(Integer, nullable=True)
    geneid = Column(Integer, nullable=True)
    sym = Column(String(20), nullable=True)
    family = Column(String(255), nullable=True)
    chr = Column(String(255), nullable=True)
    seq = Column(Text, nullable=True)
    dtoid = Column(String(13), nullable=True)
    stringid = Column(String(15), nullable=True)
    dtoclass = Column(String(255), nullable=True)
    preferred_symbol = Column(String(20), nullable=False)
    novelty = Column(DECIMAL(34, 16), nullable=True)

    __table_args__ = (
        Index("protein_idx1", "uniprot"),
        Index("protein_name_idx", "name"),
        Index("protein_sym_idx", "sym"),
        Index(
            "protein_text1_idx",
            "name", "description",
            mysql_prefix="FULLTEXT"
        ),
        Index(
            "protein_text2_idx",
            "uniprot",
            "sym",
            "stringid",
            mysql_prefix="FULLTEXT"
        )
    )

class Target(Base):
    __tablename__ = 'target'
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=False)
    ifx_id = Column(String(18), nullable=False)
    name = Column(String(255), nullable=False)
    ttype = Column(String(255), nullable=False)
    description = Column(Text)
    comment = Column(Text)
    tdl = Column(Enum(TDL))
    idg = Column(SmallInteger, nullable=False, default=0)
    fam = Column(String(255))
    famext = Column(String(255))

class T2TC(Base):
    __tablename__ = 't2tc'
    target_id = Column(Integer, ForeignKey('target.id'), primary_key=True, nullable=False)
    protein_id = Column(Integer, ForeignKey('protein.id'), primary_key=True)

class GoA(Base):
    __tablename__ = 'goa'
    protein_id = Column(Integer, ForeignKey('protein.id'), primary_key=True)
    go_id = Column(String(255), primary_key=True)
    go_term = Column(Text)
    evidence = Column(String(3), primary_key=True)
    goeco = Column(String(255), nullable=False)
    assigned_by = Column(String(50), primary_key=True)
    go_type = Column(Enum(GoType))
    go_term_text = Column(Text)

    __table_args__ = (
        Index(
            "goa_text_idx",
            "go_term",
            mysql_prefix="FULLTEXT"
        ),
    )

class TDL_info(Base):
    __tablename__ = "tdl_info"
    id = Column(Integer, primary_key=True)
    itype = Column(String(255), nullable=False)
    target_id = Column(Integer, ForeignKey('target.id'))
    protein_id = Column(Integer, ForeignKey('protein.id'))
    string_value = Column(Text)
    number_value = Column(Float)
    integer_value = Column(Integer)
    date_value = Column(Date)
    boolean_value = Column(SmallInteger)


class GeneRif(Base):
    __tablename__ = 'generif'
    id = Column(String(64), primary_key=True, unique=True)
    protein_id = Column(Integer, ForeignKey('protein.id'), nullable=False)
    gene_id = Column(Integer, nullable=False)
    text = Column(Text, nullable=False)
    date = Column(Date)

    __table_args__ = (
        Index(
            "generif_text_idx",
            "text",
            mysql_prefix="FULLTEXT"
        ),
    )

class GeneRif2Pubmed(Base):
    __tablename__ = 'generif2pubmed'
    generif_id = Column(String(64), ForeignKey('generif.id'), primary_key=True)
    pubmed_id = Column(Integer, primary_key=True)

class Ligand(Base):
    __tablename__ = 'ncats_ligands'
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=False)
    identifier = Column(String(255), nullable=False, unique=True)
    name = Column(Text, nullable=False)
    isDrug = Column(SmallInteger)
    smiles = Column(Text)
    PubChem = Column(String(255))
    ChEMBL = Column(String(255))
    guide_to_pharmacology = Column("Guide to Pharmacology", String(255))
    DrugCentral = Column(String(255))
    description = Column(Text)
    actCnt = Column(Integer)
    targetCount = Column(Integer)
    unii = Column(String(10))
    pt = Column(String(128))

    __table_args__ = (
        Index("ncats_ligands_identifier_index", "identifier"),
        Index(
            "text_search",
            "name","ChEMBL","PubChem","Guide to Pharmacology","DrugCentral",
            mysql_prefix="FULLTEXT"
        ),
    )

class LigandActivity(Base):
    __tablename__ = 'ncats_ligand_activity'
    id = Column(Integer, primary_key=True)
    ncats_ligand_id = Column(Integer, ForeignKey('ncats_ligands.id'), nullable=False)
    target_id = Column(Integer, ForeignKey('target.id'), nullable=False)
    act_value = Column(Float)
    act_type = Column(String(255))
    action_type = Column(String(255))
    reference = Column(Text)
    reference_source = Column(String(255))
    pubmed_ids = Column(Text)  # pipe delimited list]]
    __table_args__ = (
        Index("target_ligand_index", "target_id", "ncats_ligand_id"),
    )

class DrugActivity(Base):
    __tablename__ = 'drug_activity'
    id = Column(Integer, primary_key=True)
    target_id = Column(Integer, ForeignKey('target.id'), nullable=False)
    drug = Column(String(255), nullable=False)
    act_value = Column(Float)
    act_type = Column(String(255))
    action_type = Column(String(255))
    has_moa = Column(SmallInteger)
    source = Column(String(255))
    reference = Column(Text)
    smiles = Column(Text)
    cmpd_chemblid = Column(String(255))
    nlm_drug_info = Column(Text)
    cmpd_pubchem_cid = Column(Integer)
    dcid = Column(Integer)
    lychi_h4 = Column(String(15))

class PPI(Base):
    __tablename__ = "ncats_ppi"
    id = Column(Integer, primary_key=True)
    ppitypes = Column(String(255), nullable=False)
    protein_id = Column(Integer, ForeignKey('protein.id'), nullable=False)
    other_id = Column(Integer, ForeignKey('protein.id'), nullable=False)
    p_int = Column(Float)
    p_ni = Column(Float)
    p_wrong = Column(Float)
    evidence = Column(String(255))
    interaction_type = Column(String(100))
    score = Column(Integer)

class GO(Base):
    __tablename__ = 'go'

    go_id = Column(String(12), primary_key=True)  # varchar(12)
    name = Column(Text, nullable=False)
    namespace = Column(Text, nullable=False)
    def_ = Column('def', Text)  # "def" is a reserved keyword in Python
    comment = Column(Text)


class GOParent(Base):
    __tablename__ = 'go_parent'

    go_id = Column(String(12), ForeignKey('go.go_id'), nullable=False, primary_key=True)
    parent_id = Column(String(12), ForeignKey('go.go_id'), nullable=False, primary_key=True)

class MondoXref(Base):
    __tablename__ = 'mondo_xref'
    __table_args__ = (
        Index("mondo_xref_idx1", "mondoid"),
        Index("mondo_xref_id_map", "equiv_to", "db", "value"),
        Index("mondo_xref_equiv_to_xref_index", "equiv_to", "xref"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    mondoid = Column(String(20), ForeignKey("mondo.mondoid"), nullable=False)
    db = Column(String(24), nullable=False)
    value = Column(String(255), nullable=False)
    equiv_to = Column(Boolean, nullable=False, default=False)  # tinyint(1) in MySQL is usually Boolean
    source_info = Column(Text)
    xref = Column(String(255), nullable=True)

class MondoParent(Base):
    __tablename__ = "mondo_parent"

    mondoid = Column(String(20),
                     ForeignKey("mondo.mondoid"),
                     nullable=False, primary_key=True)
    parentid = Column(String(20), nullable=False, primary_key=True)

    __table_args__ = (
        Index("mondo_parent_idx1", "mondoid"),
    )

class Mondo(Base):
    __tablename__ = "mondo"
    __table_args__ = (
        Index("mondo_name_index", "name", mysql_length={"name": 191}),
    )

    mondoid = Column(String(20), primary_key=True, nullable=False)
    name = Column(Text, nullable=False)
    def_ = Column("def", Text)  # 'def' is a Python keyword
    comment = Column(Text)

class Disease(Base):
    __tablename__ = "disease"
    __table_args__ = (
        Index("disease_idx1", "dtype"),
        Index("disease_idx2", "protein_id"),
        Index("disease_idx3", "nhprotein_id"),
        Index("fk_disease_mondo", "mondoid"),
        Index("disease_idx5", "did"),
        Index("disease_idx6", "drug_name", mysql_length={"drug_name": 256}),
        Index("disease_idx4", "ncats_name", mysql_length={"ncats_name": 256}),
        Index("disease_text_idx", "ncats_name", "description", "drug_name", mysql_prefix="FULLTEXT")
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    dtype = Column(String(255),
                   ForeignKey("disease_type.name"), nullable=False)
    protein_id = Column(Integer, ForeignKey('protein.id'))
    nhprotein_id = Column(Integer, ForeignKey("nhprotein.id"))
    name = Column(Text, nullable=False)
    ncats_name = Column(Text, nullable=False)
    did = Column(String(20))
    evidence = Column(Text)
    zscore = Column(DECIMAL(5, 3))
    conf = Column(DECIMAL(2, 1))
    description = Column(Text)
    reference = Column(String(255))
    drug_name = Column(Text)
    log2foldchange = Column(DECIMAL(5, 3))
    pvalue = Column(String(255))
    score = Column(DECIMAL(16, 15))
    source = Column(String(255))
    O2S = Column(DECIMAL(16, 13))
    S2O = Column(DECIMAL(16, 13))
    mondoid = Column(String(20),
                     ForeignKey("mondo.mondoid"))
    updated = Column(Boolean)

class DiseaseType(Base):
    __tablename__ = "disease_type"

    name = Column(String(255), primary_key=True, nullable=False)
    description = Column(Text)

class NhProtein(Base):
    __tablename__ = "nhprotein"
    __table_args__ = (
        Index("nhprotein_idx1", "uniprot"),
        Index("nhprotein_idx2", "sym"),
        Index("nhprotein_idx3", "taxid", "geneid"),
        Index("nhprotein_idx4", "species")
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot = Column(String(20), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text)
    sym = Column(String(30))
    species = Column(String(40), nullable=False)
    taxid = Column(Integer, nullable=False)
    geneid = Column(Integer)

class NcatsDisease(Base):
    __tablename__ = "ncats_disease"
    __table_args__ = (
        Index("ncats_disease_mondoid_foreign", "mondoid"),
        Index("ncats_disease_gard_rare_index", "gard_rare"),
        Index("ncats_disease_name_index", "name"),
        Index(
            "ncats_disease_fulltext_idx",
            "name", "uniprot_description", "do_description", "mondo_description",
            mysql_prefix="FULLTEXT"
        )
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    ifx_id = Column(String(64), nullable=True, unique=True)
    name = Column(String(255), nullable=False)
    uniprot_description = Column(Text)
    do_description = Column(Text)
    mondo_description = Column(Text)
    mondoid = Column(String(20), ForeignKey("mondo.mondoid"))
    maxTDL = Column(String(6))
    target_count = Column(Integer)
    direct_target_count = Column(Integer)
    gard_rare = Column(Boolean)
    novelty = Column(DECIMAL(34, 16), nullable=True)

class NcatsDataSourceMap(Base):
    __tablename__ = "ncats_dataSource_map"
    __table_args__ = (
        Index("dataSource_dataSource", "dataSource"),
        Index("dataSource_protein", "protein_id"),
        Index("dataSource_ligand", "ncats_ligand_id"),
        Index("dataSource_disease", "disease_name", mysql_length={"disease_name": 200})
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataSource = Column(String(50), nullable=False)
    url = Column(String(128))
    license = Column(String(128))
    licenseURL = Column(String(128))
    protein_id = Column(Integer, ForeignKey("protein.id"))
    ncats_ligand_id = Column(Integer, ForeignKey("ncats_ligands.id"))
    disease_name = Column(Text)

class Ortholog(Base):
    __tablename__ = "ortholog"
    __table_args__ = (
        Index("ortholog_idx1", "protein_id"),
        Index("ortholog_idx2", "taxid", "geneid"),
        Index("ortholog_idx3", "symbol"),
        Index("ortholog_facet_idx", "protein_id", "geneid", "taxid"),
        Index("ortholog_text_idx", "symbol", "name", mysql_prefix="FULLTEXT")
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    taxid = Column(Integer, nullable=False)
    species = Column(String(255), nullable=False)
    db_id = Column(String(255))
    geneid = Column(Integer)
    symbol = Column(String(255), nullable=False)
    name = Column(String(255), nullable=False)
    mod_url = Column(Text)
    sources = Column(String(255), nullable=False)

class Pathway(Base):
    __tablename__ = "pathway"
    __table_args__ = (
        Index("pathway_idx1", "target_id"),
        Index("pathway_idx2", "protein_id"),
        Index("pathway_idx3", "pwtype"),
        Index("pathway_idx4", "name", mysql_length={"name": 256}),
        Index("pathway_idx5", "protein_id", "pwtype", "name", mysql_length={"name": 256}),
        Index("pathway_text_idx", "name", mysql_prefix="FULLTEXT")
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    target_id = Column(Integer, ForeignKey("target.id", ondelete="CASCADE"))
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"))
    pwtype = Column(String(255), nullable=False)
    id_in_source = Column(String(255))
    name = Column(Text, nullable=False)
    url = Column(Text)

class Phenotype(Base):
    __tablename__ = "phenotype"
    __table_args__ = (
        Index("phenotype_idx1", "ptype"),
        Index("phenotype_idx2", "protein_id"),
        Index("phenotype_idx3", "nhprotein_id"),
        Index("phenotype_idx4", "term_name"),
        Index("phenotype_idx5", "ptype", "nhprotein_id"),
        Index("phenotype_idx6", "ptype", "nhprotein_id", "protein_id"),
        Index("phenotype_nhid_idx", "nhprotein_id"),
        Index("phenotype_text_idx", "term_name", "term_description", mysql_prefix="FULLTEXT")
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    ptype = Column(String(255), nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"))
    nhprotein_id = Column(Integer, ForeignKey("nhprotein.id", ondelete="CASCADE"))
    trait = Column(Text)
    top_level_term_id = Column(String(255))
    top_level_term_name = Column(String(255))
    term_id = Column(String(255))
    term_name = Column(String(255))
    term_description = Column(Text)
    p_value = Column(Float)
    percentage_change = Column(String(255))
    effect_size = Column(String(255))
    procedure_name = Column(String(255))
    parameter_name = Column(String(255))
    gp_assoc = Column(Boolean)
    statistical_method = Column(Text)
    sex = Column(String(8))
    original_id = Column(Integer)

class Expression(Base):
    __tablename__ = "expression"
    __table_args__ = (
        Index("expression_idx1", "etype"),
        Index("expression_idx3", "protein_id"),
        Index("expression_idx4", "tissue_id"),
        Index("expression_facet_idx", "expressed", "tissue_id", "protein_id"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    etype = Column(String(255), nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"))
    source_id = Column(String(20), nullable=False)
    tissue = Column(String(255), nullable=False)
    tissue_id = Column(Integer, ForeignKey("tissue.id", ondelete="CASCADE"), nullable=False)
    qual_value = Column(Enum("Not detected", "Low", "Medium", "High"))
    number_value = Column(DECIMAL(12, 6))
    expressed = Column(Integer, nullable=False)  # tinyint(4) maps to Integer
    source_rank = Column(DECIMAL(12, 6))
    evidence = Column(String(255))
    oid = Column(String(20))
    uberon_id = Column(String(20))

class Tissue(Base):
    __tablename__ = "tissue"
    __table_args__ = (
        Index("expression_text_idx", "name", mysql_prefix="FULLTEXT"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False)

class ViralProtein(Base):
    __tablename__ = "viral_protein"

    id = Column(Integer, primary_key=True, nullable=False)
    name = Column(String(128), nullable=True)
    ncbi = Column(String(128), nullable=True)
    virus_id = Column(String(16), ForeignKey("virus.virusTaxid"), nullable=True)

    __table_args__ = (
        Index("virus_id_idx", "virus_id"),
    )

class ViralPPI(Base):
    __tablename__ = "viral_ppi"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    viral_protein_id = Column(Integer, ForeignKey("viral_protein.id"), nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id"), nullable=True)
    dataSource = Column(String(20), nullable=True)
    finalLR = Column(DECIMAL(20, 12), nullable=False)
    pdbIDs = Column(String(128), nullable=True)
    highConfidence = Column(SmallInteger, nullable=True)  # tinyint

    __table_args__ = (
        Index("viral_protein_id_idx", "viral_protein_id"),
        Index("protein_id_idx", "protein_id"),
        Index("high_conf_idx", "highConfidence"),
    )

class Virus(Base):
    __tablename__ = "virus"

    virusTaxid = Column(String(16), primary_key=True, nullable=False)
    nucleic1 = Column(String(128), nullable=True)
    nucleic2 = Column(String(128), nullable=True)
    order = Column(String(128), nullable=True)
    family = Column(String(128), nullable=True)
    subfamily = Column(String(128), nullable=True)
    genus = Column(String(128), nullable=True)
    species = Column(String(128), nullable=True)
    name = Column(String(128), nullable=True)

class Xref(Base):
    __tablename__ = "xref"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    xtype = Column(String(255), nullable=False)
    target_id = Column(Integer, ForeignKey("target.id", ondelete="CASCADE"), nullable=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=True)
    nucleic_acid_id = Column(Integer, nullable=True)  # No foreign key specified here
    value = Column(String(255), nullable=False)
    xtra = Column(String(255), nullable=True)
    dataset_id = Column(Integer, nullable=True)
    nhprotein_id = Column(Integer, ForeignKey("nhprotein.id"), nullable=True)

    __table_args__ = (
        Index("xref_idx1", "xtype"),
        Index("xref_idx2", "target_id"),
        Index("xref_idx4", "protein_id"),
        Index("xref_idx6", "dataset_id"),
        Index("fk_xref_nhprotein", "nhprotein_id"),
        Index("xref_idx7", "xtra"),
        Index("xref_type_val_idx", "xtype", "value"),
        # Note: SQLAlchemy does not have built-in support for FULLTEXT indexes.
        # You would add that manually via DDL or migrations if needed.
    )

class NihList(Base):
    __tablename__ = "nih_list"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id"), nullable=False)
    symbol = Column(String(20), nullable=False)
    family = Column(String(20), nullable=False)
    nih_list = Column(Integer, ForeignKey("nih_list_type.id"), nullable=False)

    __table_args__ = (
        Index("nih_list_protein_id_foreign", "protein_id"),
        Index("nih_list_nih_list_foreign", "nih_list"),
    )

class PantherClass(Base):
    __tablename__ = "panther_class"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    pcid = Column(String(7), nullable=False, unique=True)
    parent_pcids = Column(String(255), nullable=True)
    name = Column(Text, nullable=False)
    description = Column(Text, nullable=True)

    __table_args__ = (
        UniqueConstraint("pcid", name="panther_class_idx1"),
    )


class P2PC(Base):
    __tablename__ = "p2pc"

    panther_class_id = Column(Integer, ForeignKey("panther_class.id"), primary_key=True, nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), primary_key=True, nullable=False)

    __table_args__ = (
        Index("p2pc_idx1", "panther_class_id"),
        Index("p2pc_idx2", "protein_id"),
    )


class NihListType(Base):
    __tablename__ = "nih_list_type"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    list_type = Column(String(255), nullable=False, unique=True)

class DTO(Base):
    __tablename__ = "dto"

    dtoid = Column(String(255), primary_key=True, nullable=False)
    name = Column(Text, nullable=False)
    parent_id = Column(String(255), nullable=True)
    def_ = Column("def", Text, nullable=True)  # 'def' is a Python keyword, so renamed to def_

    __table_args__ = (
        Index("dto_idx1", "parent_id"),
    )


class DTOParent(Base):
    __tablename__ = "dto_parent"

    dtoid = Column(String(255), ForeignKey("dto.dtoid"), primary_key=True, nullable=False)
    parent_id = Column(String(255), ForeignKey("dto.dtoid"), primary_key=True, nullable=False)

    __table_args__ = (
        Index("dto_parent_idx1", "dtoid"),
        Index("dto_parent_idx2", "parent_id"),
    )


class P2DTO(Base):
    __tablename__ = "p2dto"

    dtoid = Column(String(255), ForeignKey("dto.dtoid"), primary_key=True, nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id"), primary_key=True, nullable=False)

    __table_args__ = (
        Index("p2dto_dtoid_foreign", "dtoid"),
        Index("p2dto_protein_id_foreign", "protein_id"),
    )

class NcatsD2DA(Base):
    __tablename__ = "ncats_d2da"

    ncats_disease_id = Column(Integer, ForeignKey("ncats_disease.id"), primary_key=True, nullable=False)
    disease_assoc_id = Column(Integer, ForeignKey("disease.id"), primary_key=True, nullable=False)
    direct = Column(SmallInteger, nullable=True)  # tinyint(1) as SmallInteger

    __table_args__ = (
        Index("ncats_d2da_ncats_disease_id_foreign", "ncats_disease_id"),
        Index("ncats_d2da_disease_assoc_id_foreign", "disease_assoc_id"),
    )

class Tiga(Base):
    __tablename__ = "tiga"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    ensg = Column(String(15), nullable=False)
    efoid = Column(String(15), nullable=False)
    trait = Column(String(255), nullable=False)
    n_study = Column(Integer, nullable=True)
    n_snp = Column(Integer, nullable=True)
    n_snpw = Column(DECIMAL(6, 3), nullable=True)
    geneNtrait = Column(Integer, nullable=True)
    geneNstudy = Column(Integer, nullable=True)
    traitNgene = Column(Integer, nullable=True)
    traitNstudy = Column(Integer, nullable=True)
    pvalue_mlog_median = Column(DECIMAL(7, 3), nullable=True)
    pvalue_mlog_max = Column(DECIMAL(8, 3), nullable=True)
    or_median = Column(DECIMAL(8, 3), nullable=True)
    n_beta = Column(Integer, nullable=True)
    study_N_mean = Column(Integer, nullable=True)
    rcras = Column(DECIMAL(5, 3), nullable=True)
    meanRank = Column(DECIMAL(18, 12), nullable=True)
    meanRankScore = Column(DECIMAL(18, 14), nullable=True)
    ncats_disease_id = Column(Integer, ForeignKey("ncats_disease.id"), nullable=True)

    __table_args__ = (
        Index("tiga_idx1", "protein_id"),
        Index("tiga_trait_protein", "trait", "protein_id"),
        Index("tiga_ncats_disease_id_foreign", "ncats_disease_id"),
    )


class TigaProvenance(Base):
    __tablename__ = "tiga_provenance"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    ensg = Column(String(15), nullable=False)
    efoid = Column(String(15), nullable=False)
    study_acc = Column(String(20), nullable=False)
    pubmedid = Column(Integer, nullable=False)

class DrgcResource(Base):
    __tablename__ = "drgc_resource"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    rssid = Column(Text, nullable=False)
    resource_type = Column(String(255), nullable=False)
    target_id = Column(Integer, ForeignKey("target.id", ondelete="CASCADE"), nullable=False)
    json = Column(Text, nullable=False)

    __table_args__ = (
        Index("drgc_resource_idx1", "target_id"),
    )

class TinxImportance(Base):
    __tablename__ = "tinx_importance"

    ncats_disease_id = Column(Integer, ForeignKey("ncats_disease.id"), primary_key=True, nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id"), primary_key=True, nullable=False)
    doid = Column(String(20), nullable=True)
    score = Column(DECIMAL(34, 16), nullable=False)

    __table_args__ = (
        Index("tinx_importance_idx1", "protein_id"),
        Index("tinx_importance_idx2", "ncats_disease_id"),
        Index("tinx_importance_idx3", "doid"),
    )

class GeneAttributeType(Base):
    __tablename__ = "gene_attribute_type"

    id = Column(Integer, primary_key=True, autoincrement=True, nullable=False)
    name = Column(String(255), nullable=False, unique=True)
    association = Column(Text, nullable=False)
    description = Column(Text, nullable=False)
    resource_group = Column(
        Enum(
            'omics',
            'genomics',
            'proteomics',
            'physical interactions',
            'transcriptomics',
            'structural or functional annotations',
            'disease or phenotype associations',
            name="resource_group_enum"
        ),
        nullable=False
    )
    measurement = Column(String(255), nullable=False)
    attribute_group = Column(String(255), nullable=False)
    attribute_type = Column(String(255), nullable=False)
    pubmed_ids = Column(Text, nullable=True)
    url = Column(Text, nullable=True)

class DO(Base):
    __tablename__ = "do"

    doid = Column(String(20), primary_key=True, nullable=False)
    name = Column(Text, nullable=False)
    def_ = Column("def", Text, nullable=True)  # 'def' is a reserved keyword in Python


class DOParent(Base):
    __tablename__ = "do_parent"

    doid = Column(String(20), ForeignKey("do.doid"), nullable=False, primary_key=True)
    parent_id = Column(String(20), nullable=False, primary_key=True)

    __table_args__ = (
        Index("fk_do_parent__do", "doid"),
    )

class Protein2Pubmed(Base):
    __tablename__ = "protein2pubmed"

    id = Column(Integer, primary_key=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    pubmed_id = Column(Integer, nullable=False)
    gene_id = Column(Integer, nullable=True)
    source = Column(String(45), nullable=False)

    __table_args__ = (
        Index("protein2pubmed_idx1", "protein_id"),
        Index("protein2pubmed_idx2", "pubmed_id"),
        Index("protein2pubmed_type", "source"),
        Index("protein2pubmed_type_protein_id", "source", "protein_id"),
        Index("protein2pubmed_unique_assoc", "protein_id", "pubmed_id", "gene_id", "source", unique=True),
    )

class Alias(Base):
    __tablename__ = "alias"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    type = Column(ENUM('symbol','uniprot','NCBI Gene ID'), nullable=False)
    value = Column(String(255), nullable=False)
    dataset_id = Column(Integer, nullable=True)

    __table_args__ = (
        Index("alias_idx1", "protein_id"),
        Index("alias_idx2", "dataset_id"),
        Index("alias_text_idx", "value", mysql_prefix="FULLTEXT"),
    )

class Ptscore(Base):
    __tablename__ = "ptscore"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    year = Column(Integer, nullable=False)
    score = Column(DECIMAL(12, 6), nullable=False)

    __table_args__ = (
        Index("ptscore_idx1", "protein_id"),
        Index("ptscore_score_idx", "score"),
    )

class Pmscore(Base):
    __tablename__ = "pmscore"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    year = Column(Integer, nullable=False)
    score = Column(DECIMAL(12, 6), nullable=False)

    __table_args__ = (
        Index("pmscore_idx1", "protein_id"),
    )

class PatentCount(Base):
    __tablename__ = "patent_count"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    year = Column(Integer, nullable=False)
    count = Column(Integer, nullable=False)

    __table_args__ = (
        Index("patent_count_idx1", "protein_id"),
    )

class Omim(Base):
    __tablename__ = "omim"

    mim = Column(Integer, primary_key=True)
    title = Column(String(255), nullable=False)

class InputVersion(Base):
    __tablename__ = 'input_version'

    source_key = Column(String(45), nullable=False, primary_key=True)
    data_source = Column(String(45), nullable=False)
    file_key = Column(String(45), nullable=False, primary_key=True)
    file = Column(String(256), nullable=False)
    version = Column(String(45))
    release_date = Column(Date)
    download_date = Column(Date, nullable=False)


class DataSourceVersion(Base):
    __tablename__ = "data_source_version"
    __table_args__ = (
        Index("data_source_version_source_idx", "data_source"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    data_source = Column(String(255), nullable=False)
    version = Column(String(255))
    version_date = Column(Date)
    download_date = Column(Date)
    collections = Column(Text)


class ETLRun(Base):
    __tablename__ = "etl_run"
    __table_args__ = (
        Index("etl_run_stage_idx", "stage"),
        Index("etl_run_database_idx", "database_name"),
        Index("etl_run_run_date_idx", "run_date"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    stage = Column(String(64), nullable=False)
    database_name = Column(String(255), nullable=False)
    run_date = Column(DateTime, nullable=False)
    git_commit = Column(String(255))
    git_branch = Column(String(255))
    git_tag = Column(String(255))
    runner = Column(String(255))
    hostname = Column(String(255))
    platform = Column(String(255))
    platform_version = Column(String(255))
    python_version = Column(String(255))


class ETLAdapterRun(Base):
    __tablename__ = "etl_adapter_run"
    __table_args__ = (
        Index("etl_adapter_run_run_idx", "run_id"),
        Index("etl_adapter_run_database_idx", "database_name"),
        Index("etl_adapter_run_status_idx", "status"),
        Index("etl_adapter_run_position_idx", "run_id", "adapter_position"),
        UniqueConstraint("run_id", "adapter_name", name="etl_adapter_run_unique"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(String(255), nullable=False)
    database_name = Column(String(255), nullable=False)
    adapter_name = Column(String(255), nullable=False)
    adapter_position = Column(Integer)
    adapter_total = Column(Integer)
    status = Column(String(32), nullable=False)
    records_written = Column(Integer, nullable=False, default=0)
    store_calls = Column(Integer, nullable=False, default=0)
    input_object_count = Column(Integer, nullable=False, default=0)
    converted_object_count = Column(Integer, nullable=False, default=0)
    inserted_row_count = Column(Integer, nullable=False, default=0)
    conversion_seconds = Column(Float, nullable=False, default=0)
    serialization_seconds = Column(Float, nullable=False, default=0)
    insert_seconds = Column(Float, nullable=False, default=0)
    total_store_seconds = Column(Float, nullable=False, default=0)
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    failed_at = Column(DateTime)
    error_message = Column(Text)

class NcatsDataSource(Base):
    __tablename__ = 'ncats_dataSource'

    dataSource = Column(String(50), primary_key=True, nullable=False)
    dataSourceDescription = Column(String(512))
    url = Column(String(128))
    license = Column(String(128))
    licenseURL = Column(String(128))
    citation = Column(String(512))

class SequenceAnnotation(Base):
    __tablename__ = "sequence_annotation"
    __table_args__ = (
        Index("annotation_source_protein_idx", "dataSource", "protein_id"),
        Index("fk_annotation_protein", "protein_id"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataSource = Column(String(255), nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id"), nullable=False)
    residue_start = Column(Integer, nullable=False)
    residue_end = Column(Integer, nullable=False)
    type = Column(Enum(
        'Activation Loop','Activation Segment','alphaC-beta4 Loop',
        'CMGC Insert','Gatekeeper','Linker','KeyAA','Motif','beta-strand',
        'alpha-helix','C-Lobe','C-Spine','RD Pocket','Catalytic Loop',
        'Glycine Loop','N-Lobe','R-Spine','R-Spine Shell','Subdomain'
    ), nullable=False)
    name = Column(String(255), nullable=False)


class SequenceVariant(Base):
    __tablename__ = "sequence_variant"
    __table_args__ = (
        Index("variant_source_protein_idx", "dataSource", "protein_id"),
        Index("fk_variant_protein", "protein_id"),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    dataSource = Column(String(255), nullable=False)
    protein_id = Column(Integer, ForeignKey("protein.id"), nullable=False)
    residue = Column(Integer, nullable=False)
    variant = Column(String(1), nullable=False)
    bits = Column(Float(12, 11), nullable=False)

class AncestryDO(Base):
    __tablename__ = 'ancestry_do'

    oid = Column(String(12), ForeignKey('do.doid'), index=True, primary_key=True)
    ancestor_id = Column(String(12), ForeignKey('do.doid'), index=True, primary_key=True)

    __table_args__ = (
        Index('ancestry_do_oid_foreign', 'oid'),
        Index('ancestry_do_ancestor_id_foreign', 'ancestor_id'),
    )

class AncestryDTO(Base):
    __tablename__ = 'ancestry_dto'

    oid = Column(String(255), ForeignKey('dto.dtoid'), index=True, primary_key=True)
    ancestor_id = Column(String(255), ForeignKey('dto.dtoid'), index=True, primary_key=True)

    __table_args__ = (
        Index('ancestry_dto_oid_foreign', 'oid'),
        Index('ancestry_dto_ancestor_id_foreign', 'ancestor_id'),
    )

class AncestryMONDO(Base):
    __tablename__ = 'ancestry_mondo'

    oid = Column(String(20), ForeignKey('mondo.mondoid'), index=True, primary_key=True)
    ancestor_id = Column(String(20), ForeignKey('mondo.mondoid'), index=True, primary_key=True)

    __table_args__ = (
        Index('ancestry_mondo_oid_foreign', 'oid'),
        Index('ancestry_mondo_ancestor_id_foreign', 'ancestor_id'),
    )

class AncestryUBERON(Base):
    __tablename__ = 'ancestry_uberon'

    oid = Column(String(30), ForeignKey('uberon.uid'), index=True, primary_key=True)
    ancestor_id = Column(String(30), ForeignKey('uberon.uid'), index=True, primary_key=True)

    __table_args__ = (
        Index('ancestry_uberon_oid_foreign', 'oid'),
        Index('ancestry_uberon_ancestor_id_foreign', 'ancestor_id'),
    )

class Uberon(Base):
    __tablename__ = "uberon"

    uid = Column(String(30), primary_key=True)
    name = Column(Text, nullable=False)
    def_ = Column("def", Text)
    comment = Column(Text)

    __table_args__ = (
        Index(
            "uberon_text_idx",
            "name",
            "def",
            "comment",
            mysql_prefix="FULLTEXT"
        ),
    )


class UberonParent(Base):
    __tablename__ = "uberon_parent"

    uid = Column(String(30),
                 ForeignKey("uberon.uid"), nullable=False, index=True, primary_key=True)
    parent_id = Column(String(30), nullable=False, primary_key=True)


class UberonXref(Base):
    __tablename__ = "uberon_xref"

    uid = Column(String(30),
                 ForeignKey("uberon.uid"), nullable=False, primary_key=True)
    db = Column(String(24), nullable=False, primary_key=True)
    value = Column(String(255), nullable=False, primary_key=True)

class KeggNearestTclin(Base):
    __tablename__ = "kegg_nearest_tclin"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False, index=True)
    tclin_id = Column(Integer, ForeignKey("target.id", ondelete="CASCADE"), nullable=False, index=True)
    direction = Column(Enum('upstream', 'downstream'), nullable=False)
    distance = Column(Integer, nullable=False)

    __table_args__ = (
        Index("kegg_nearest_tclin_idx1", "protein_id"),
        Index("kegg_nearest_tclin_idx2", "tclin_id"),
    )

class KeggDistance(Base):
    __tablename__ = "kegg_distance"

    id = Column(Integer, primary_key=True, autoincrement=True)
    pid1 = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False, index=True)
    pid2 = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False, index=True)
    distance = Column(Integer, nullable=False)

    __table_args__ = (
        Index("kegg_distance_idx1", "pid1"),
        Index("kegg_distance_idx2", "pid2"),
    )

class Gtex(Base):
    __tablename__ = "gtex"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), index=True, nullable=True)
    tissue = Column(Text, nullable=False)
    tpm = Column(DECIMAL(12, 6), nullable=False)
    tpm_rank = Column(DECIMAL(4, 3), nullable=True)
    tpm_male = Column(DECIMAL(12, 6), nullable=True)
    tpm_male_rank = Column(DECIMAL(4, 3), nullable=True)
    tpm_female = Column(DECIMAL(12, 6), nullable=True)
    tpm_female_rank = Column(DECIMAL(4, 3), nullable=True)
    uberon_id = Column(String(30), ForeignKey("uberon.uid"), index=True, nullable=True)

    __table_args__ = (
        Index("expression_idx1", "protein_id"),
        Index("fk_gtex_uberon", "uberon_id"),
    )


class ExtLink(Base):
    __tablename__ = "extlink"

    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(Integer, ForeignKey("protein.id", ondelete="CASCADE"), nullable=False)
    source = Column(String(255), ForeignKey("affiliate.source"), nullable=False)
    url = Column(Text, nullable=False)

    # Index
    __table_args__ = (
        Index("extlink_idx1", "protein_id"),
    )

class Affiliate(Base):
    __tablename__ = "affiliate"

    id = Column(Integer, primary_key=True, autoincrement=True)
    source = Column(String(255), nullable=False, unique=True)
    display_name = Column(String(255), nullable=False)
    description = Column(String(255), nullable=False)
    link_count = Column(Integer, nullable=False, default=0)

class NcatsTypeaheadIndex(Base):
    __tablename__ = "ncats_typeahead_index"

    id = Column(Integer, primary_key=True, autoincrement=True)
    category = Column(String(45), nullable=False)
    value = Column(String(255), nullable=False)
    reference_id = Column(String(255))

    __table_args__ = (
        Index("ncats_typeahead_text", "value", mysql_prefix="FULLTEXT"),
    )


class WordCount(Base):
    __tablename__ = "word_count"

    word = Column(String(128), primary_key=True)
    count = Column(Integer)
