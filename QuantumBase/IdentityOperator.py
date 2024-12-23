from .QuantumState import QuantumState
from .Operator import Operator

class IdentityOperator(Operator):
    def __init__(self):
        super().__init__()