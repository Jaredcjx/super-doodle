import sympy as smp
from .QuantumState import QuantumState
from sympy.physics.quantum.constants import hbar

class MomentumQuantumState(QuantumState):
    #Constructor
    def __init__(self, function, lowerLimit=None, upperLimit=None):
        super().__init__(function, lowerLimit, upperLimit)
    
    #Methods
    def __mul__(self, other) -> complex:
        if (not isinstance(other, MomentumQuantumState)) and (isinstance(other, QuantumState)):
            raise TypeError('Convert argument to Momentum Basis First')
        else:
            return super().__mul__(other)
    
    def update_state(self, wavefunction : smp.Function):
        return MomentumQuantumState(wavefunction, self.lowerLimit, self.upperLimit)
    
    def get_position_basis_wavefunction(self):
        t, p = smp.symbols('t p')
        phi_p = lambda x: smp.exp(smp.I * p * x / hbar) / smp.sqrt(2 * smp.pi * hbar)
        position_function = smp.integrate(self.wavefunction(p) * phi_p(t), (p, -smp.oo, smp.oo), meijerg=True).simplify()
        res = lambda x: position_function.subs(t, x)
        return res