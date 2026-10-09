from typing import List, Generator
from src.constants import Prefix, DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.gene import Gene
from src.models.node import EquivalentId
from src.models.transcript import GeneTranscriptEdge, ParentLinkDetail, Transcript
from src.shared.targetgraph_parser import TargetGraphTranscriptParser


class GeneTranscriptEdgeAdapter(InputAdapter, TargetGraphTranscriptParser):
    def __init__(self, data_source):
        self.version_info = data_source.version_info()
        file_path = str(data_source.file("transcript_ids.tsv"))
        TargetGraphTranscriptParser.__init__(self, file_path=file_path)

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.TargetGraph

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[List[GeneTranscriptEdge], None, None]:
        relationships = []

        for line in self.all_rows():
            transcript_id = TargetGraphTranscriptParser.get_id(line)
            transcript_obj = Transcript(id=transcript_id)

            created = TargetGraphTranscriptParser.get_creation_date(line)
            updated = TargetGraphTranscriptParser.get_updated_time(line)

            parent_id = TargetGraphTranscriptParser.get_parent_gene_id(line)
            if 'parent_ifx_gene_id' in self.fieldnames:
                if not parent_id:
                    raise ValueError(f"Transcript {transcript_id} has no parent_ifx_gene_id")
                gene_id = parent_id
                details = [ParentLinkDetail(
                    resolution_source=line['parent_resolution_source'],
                    confidence=line['parent_confidence'],
                )]
            else:
                ensg_id = TargetGraphTranscriptParser.get_associated_ensg_id(line)
                ncbi_id = TargetGraphTranscriptParser.get_associated_ncbi_id(line)
                if ensg_id:
                    gene_id = EquivalentId(id=ensg_id, type=Prefix.ENSEMBL).id_str()
                elif ncbi_id:
                    gene_id = EquivalentId(id=ncbi_id, type=Prefix.NCBIGene).id_str()
                else:
                    raise ValueError(f"Transcript {transcript_id} has no associated gene")
                details = []

            relationships.append(
                GeneTranscriptEdge(
                    start_node=Gene(id=gene_id),
                    end_node=transcript_obj,
                    created=created,
                    updated=updated,
                    details=details,
                )
            )

        yield relationships
