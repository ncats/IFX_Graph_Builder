from typing import Generator, List

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import ReactionClass as SqliteReactionClass
from src.models.reaction import ReactionClass, ReactionClassParentEdge


def get_parent_ec(ec: str):
    root, grand, parent, leaf = ec.split('.')
    if grand == '-':
        return None
    if parent == '-':
        return f'{root}.-.-.-'
    if leaf == '-':
        return f'{root}.{grand}.-.-'
    return f'{root}.{grand}.{parent}.-'


class ReactionClassEdgeAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[ReactionClassParentEdge], None, None]:
        results = self.get_session().query(
            SqliteReactionClass.rxn_class_ec
        ).filter(SqliteReactionClass.ec_level == 4).distinct().all()

        relationship_set = set()

        for row in results:
            ec = row[0]
            parent = get_parent_ec(ec)
            while parent is not None:
                relationship_set.add(f"{ec}|{parent}")
                ec = parent
                parent = get_parent_ec(ec)

        relationships: List[ReactionClassParentEdge] = [
            ReactionClassParentEdge(
                start_node=ReactionClass(id=pair.split('|')[0]),
                end_node=ReactionClass(id=pair.split('|')[1])
            ) for pair in relationship_set
        ]
        yield relationships

