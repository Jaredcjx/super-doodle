import numpy as np

class CompQuantumState:
    #attributes
    basisOrder : int
    coefficients : np.ndarray

    #Constructor
    def __init__(self, basisOrder : int, coefficients : np.ndarray):
        self.basisOrder = basisOrder
        self.coefficients = coefficients.reshape(-1, 1)
    
    #Methods
    def left_overlap(self, other : 'CompQuantumState') -> complex:
        return np.vdot(self.coefficients, other.coefficients)
    
    def right_overlap(self, other : 'CompQuantumState') -> complex:
        return np.vdot(other.coefficients, self.coefficients)
    