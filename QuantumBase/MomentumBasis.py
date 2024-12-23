import sympy as smp
from .MomentumQuantumState import MomentumQuantumState

class MomentumBasis(MomentumQuantumState):
    #Constructor
    def __init__(self, function : smp.DiracDelta, lowerLimit=None, upperLimit=None):
        super().__init__(function, lowerLimit, upperLimit)