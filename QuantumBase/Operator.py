from .QuantumState import QuantumState
import sympy as smp

class Operator():
    operation : smp.Function # A function that takes a function as input and outputs a function

    def __init__(self, f : smp.Function = None):
        if f == None:
            self.operation = lambda f: f # Identity function by default
        else:
            self.operation = f

    def __mul__(self, other):
        if isinstance(other, Operator):
            return Operator(lambda f: self.operation(other.operation(f)))
        return Operator(lambda f: lambda x: smp.simplify(self.operation(f)(x) * other))
        
    def __add__(self, other):
        if not isinstance(other, Operator):
            raise TypeError('Can only add operators with operators')
        return Operator(lambda f: lambda x: smp.simplify(self.operation(f)(x) + other.operation(f)(x)))
    
    def __sub__(self, other):
        if not isinstance(other, Operator):
            raise TypeError('Can only subtract operators with operators')
        return Operator(lambda f: lambda x: smp.simplify(self.operation(f)(x) - other.operation(f)(x)))
    
    def __pow__(self, num : int):
        if not isinstance(num, int) and not num >= 0:
            raise TypeError('Must be positive integer')
        res = Operator()
        for _ in range(0, num):
            res = res * self
        return res
    
    def commutator(O1, O2):
        return O1 * O2 - O2 * O1
    
    def apply(self, state : QuantumState) -> QuantumState:
        wavefunction = lambda t: smp.simplify(self.operation(state.wavefunction)(t))
        return state.update_state(wavefunction)