from .QuantumState import QuantumState
from sympy.physics.quantum.constants import hbar
import sympy as smp

class PositionQuantumState(QuantumState):
    #Constructor
    def __init__(self, function, lowerLimit=None, upperLimit=None):
        super().__init__(function, lowerLimit, upperLimit)
    
    #Methods
    def __mul__(self, other) -> complex:
        if (not isinstance(other, PositionQuantumState)) and (isinstance(other, QuantumState)):
            raise TypeError('Convert argument to Position Basis First')
        else:
            return super().__mul__(other)
        
    
    def update_state(self, wavefunction : smp.Function):
        return PositionQuantumState(wavefunction, self.lowerLimit, self.upperLimit)
    
    def get_momentum_basis_wavefunction(self):
        t, x = smp.symbols('t x')
        conjugate_phi_p = lambda x: smp.exp(-smp.I * t * x / hbar) / smp.sqrt(2 * smp.pi * hbar)
        momentum_function = smp.integrate(self.wavefunction(x) * conjugate_phi_p(x), (x, -smp.oo, smp.oo), meijerg=True).simplify()
        res = lambda p: momentum_function.subs(t, p)
        return res