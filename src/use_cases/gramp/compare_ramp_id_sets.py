from typing import List
import yaml
from src.core.etl import ETL
from src.core.registry_integration import RegistryIntegration
from src.input_adapters.sqlite_ramp._compare_id_sets import MetaboliteSetEdgeAdapter
from src.input_adapters.sqlite_ramp.analyte_pathway_edge_adapter import MetabolitePathwayEdgeAdapter
from src.input_adapters.sqlite_ramp.analyte_synonym_adapter import MetaboliteSynonymAdapter
from src.input_adapters.sqlite_ramp.metabolite_adapter import MetaboliteAdapter
from src.input_adapters.sqlite_ramp.metabolite_chem_props_adapter import MetaboliteChemPropsAdapter
from src.input_adapters.sqlite_ramp.metabolite_class_adapter import MetaboliteClassAdapter
from src.input_adapters.sqlite_ramp.metabolite_class_edge_adapter import MetaboliteClassEdgeAdapter
from src.input_adapters.sqlite_ramp.pathway_adapter import PathwayAdapter
from src.interfaces.input_adapter import InputAdapter

from src.interfaces.labeler import AuxLabeler, ComparingLabeler, RaMPLabeler
from src.output_adapters.cypher_output_adapter import MemgraphOutputAdapter
from src.shared.db_credentials import DBCredentials

credentials_file = "./src/use_cases/secrets/ifxdev_memgraph.yaml"

with open(credentials_file, "r") as file:
    credentials = yaml.safe_load(file)

class build_db_for_comparing_ramp_ids:
    left_db: str
    left_label: str
    right_db: str
    right_label: str
    third_db: str
    third_label: str
    output_adapter: MemgraphOutputAdapter
    etl: ETL
    def __init__(self, left_db: str, left_label: str,
                 right_db: str, right_label: str,
                 third_db: str = None, third_label: str = None):
        self.left_db = left_db
        self.left_label = left_label
        self.right_db = right_db
        self.right_label = right_label
        self.third_db = third_db
        self.third_label = third_label

        self.output_adapter = MemgraphOutputAdapter(credentials=DBCredentials(
            url = credentials['url'],
            user = credentials['user'],
            password=credentials['password'],
            port=credentials['port'],
        ))

        self.etl = ETL(input_adapters=[], output_adapters=[self.output_adapter])


    def truncate_old_db(self):
        self.etl.create_or_truncate_datastores()

    def do_etl(self):
        # Common Input Adapters
        pathway_adapter = PathwayAdapter(data_source=self.left_db).set_single_source(True)
        metabolite_class_adapter = MetaboliteClassAdapter(data_source=self.left_db).set_single_source(True)

        self.do_etl_step([
            pathway_adapter,
            metabolite_class_adapter
        ])

        # Input Adapters
        left_met_adapter = MetaboliteAdapter(data_source=self.left_db).set_single_source(True)
        left_metabolite_synonym_list_adapter = MetaboliteSynonymAdapter(data_source=self.left_db)
        left_met_chem_props_adapter = MetaboliteChemPropsAdapter(data_source=self.left_db).set_single_source(True)
        left_metabolite_pathway_edge_adapter = MetabolitePathwayEdgeAdapter(data_source=self.left_db).set_single_source(True)
        left_metabolite_class_edge_adapter = MetaboliteClassEdgeAdapter(data_source=self.left_db).set_single_source(True)

        left_labeler = AuxLabeler(self.left_label)

        self.do_etl_step([
            left_met_adapter,
            left_metabolite_synonym_list_adapter,
            left_met_chem_props_adapter,
            left_metabolite_pathway_edge_adapter,
            left_metabolite_class_edge_adapter
        ], left_labeler)
        # raise Exception('Debug stop')
        right_met_adapter = MetaboliteAdapter(data_source=self.right_db).set_single_source(True)
        right_metabolite_synonym_list_adapter = MetaboliteSynonymAdapter(data_source=self.right_db)
        right_met_chem_props_adapter = MetaboliteChemPropsAdapter(data_source=self.right_db).set_single_source(True)
        right_metabolite_pathway_edge_adapter = MetabolitePathwayEdgeAdapter(data_source=self.right_db).set_single_source(True)
        right_metabolite_class_edge_adapter = MetaboliteClassEdgeAdapter(data_source=self.right_db).set_single_source(True)

        right_labeler = AuxLabeler(self.right_label)

        self.do_etl_step([
            right_met_adapter,
            right_metabolite_synonym_list_adapter,
            right_met_chem_props_adapter,
            right_metabolite_pathway_edge_adapter,
            right_metabolite_class_edge_adapter
        ], right_labeler)

        metabolite_set_edge_adapter = (MetaboliteSetEdgeAdapter()
                                               .set_left(self.left_db).set_right(self.right_db).set_single_source(True))
        self.do_etl_step([
            metabolite_set_edge_adapter
        ], ComparingLabeler().set_left_labeler(left_labeler).set_right_labeler(right_labeler))

        if self.third_db:
            third_met_adapter = MetaboliteAdapter(data_source=self.third_db)
            third_metabolite_synonym_list_adapter = MetaboliteSynonymAdapter(data_source=self.third_db)
            third_met_chem_props_adapter = MetaboliteChemPropsAdapter(data_source=self.third_db)
            third_metabolite_pathway_edge_adapter = MetabolitePathwayEdgeAdapter(data_source=self.third_db)
            third_metabolite_class_edge_adapter = MetaboliteClassEdgeAdapter(data_source=self.third_db)

            third_labeler = AuxLabeler(self.third_label)

            self.do_etl_step([
                third_met_adapter,
                third_metabolite_synonym_list_adapter,
                third_met_chem_props_adapter,
                third_metabolite_pathway_edge_adapter,
                third_metabolite_class_edge_adapter
            ], third_labeler)

            metabolite_set_edge_adapter = MetaboliteSetEdgeAdapter().set_left(self.left_db).set_right(self.third_db)
            self.do_etl_step([
                metabolite_set_edge_adapter,
            ], ComparingLabeler().set_left_labeler(left_labeler).set_right_labeler(third_labeler)
            )

            metabolite_set_edge_adapter = MetaboliteSetEdgeAdapter().set_left(self.right_db).set_right(self.third_db)
            self.do_etl_step([
                metabolite_set_edge_adapter
            ], ComparingLabeler().set_left_labeler(right_labeler).set_right_labeler(third_labeler)
            )


    def do_etl_step(self, input_list: List[InputAdapter], labeler: RaMPLabeler = RaMPLabeler()):
        self.etl.set_labeler(labeler)
        self.etl.input_adapters = input_list
        self.etl.do_etl()

registry = RegistryIntegration.connect({"credentials": "./src/use_cases/secrets/aws_ifx_registry.yaml"})
released_ramp = registry.resolve("ramp:sqlite_database:3.0.7")
new_ramp = registry.resolve("ramp:sqlite_database:3.0.12-ramp4.0")

build_engine = build_db_for_comparing_ramp_ids(released_ramp, "ramp3011", new_ramp, "ramp3012")
build_engine.truncate_old_db()
build_engine.do_etl()
