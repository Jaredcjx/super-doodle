import sympy as smp
from sympy.physics.quantum.constants import hbar
from .Operator import Operator

h = smp.symbols('h', real = True, positive = True)

class MomentumOperator(Operator):
    def __init__(self):
        self.operation = lambda f: lambda x: smp.simplify(-smp.I * hbar * smp.diff(f(x), x))