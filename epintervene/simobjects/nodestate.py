from enum import Enum, unique, auto

@unique
class NodeState(Enum):
    SUSCEPTIBLE = auto()
    EXPOSED = auto()
    INFECTED = auto()
    RECOVERED = auto()
    VACCINATED = auto()