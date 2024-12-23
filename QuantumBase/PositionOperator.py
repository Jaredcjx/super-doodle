from .Operator import Operator
import sympy as smp

class PositionOperator(Operator):
    def __init__(self):
        self.operation = lambda f: lambda x: smp.simplify(x * f(x))