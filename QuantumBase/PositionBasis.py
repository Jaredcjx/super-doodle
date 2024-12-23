import sympy as smp
from .PositionQuantumState import PositionQuantumState

class PositionBasis(PositionQuantumState):
    #Constructor
    def __init__(self, function : smp.DiracDelta, lowerLimit=None, upperLimit=None):
        super().__init__(function, lowerLimit, upperLimit)