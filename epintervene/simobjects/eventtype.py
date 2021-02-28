from enum import Enum, unique, auto

@unique
class EventType(Enum):
    INFECTEDSUSCEPTIBLE = auto()
    EXPOSEDSUSCEPTIBLE = auto()
    RECOVER = auto()
    EXPOSEDINFECTED = auto()