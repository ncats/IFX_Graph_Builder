from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

from src.core.decorators import facets
from src.models.node import Node, Relationship
from src.models.pounce.demographics import Demographics

if TYPE_CHECKING:
    from src.models.pounce.biospecimen import Biospecimen


@dataclass
@facets(category_fields=['type'])
class Biosample(Node):
    original_id: str = None
    type: str = None
    tracking: Optional[str] = None
    demographics: Optional[Demographics] = None
    geo_characteristics: Optional[str] = None


@dataclass
class BiosampleBiospecimenEdge(Relationship):
    start_node: "Biosample" = None
    end_node: "Biospecimen" = None
