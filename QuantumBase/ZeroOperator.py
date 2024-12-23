from .QuantumState import QuantumState
from .Operator import Operator

class ZeroOperator(Operator):
    def apply(self, state : QuantumState) -> QuantumState:
        return state.update_state(lambda x: 0)