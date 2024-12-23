import sympy as smp
class QuantumState:

    wavefunction : smp.Function
    lowerLimit = -smp.oo
    upperLimit = smp.oo

    def __init__(self, function : smp.Function, lowerLimit = None, upperLimit = None):
        if lowerLimit is not None and upperLimit is not None:
            self.lowerLimit = lowerLimit
            self.upperLimit = upperLimit
        self.wavefunction = function

    #Methods
    def normalised_wavefunction(self, numerical : bool = None):
        if numerical == None:
            normalised_wf = QuantumState.__normalise_wavefunction_exact(self.wavefunction, self.lowerLimit, self.upperLimit)
        elif numerical:
            normalised_wf = QuantumState.__normalise_wavefunction_numer(self.wavefunction, self.lowerLimit, self.upperLimit)
        return normalised_wf

    def __mul__(self, other):
        if isinstance(other, QuantumState): # <self|other>
            def helper(t : smp.Symbol):
                normalised_wf1 = self.normalised_wavefunction()
                print('Left normalised wavefunction is', normalised_wf1(t))
                normalised_wf2 = other.normalised_wavefunction()
                print('Right normalised wavefunction is', normalised_wf2(t))
                return smp.integrate(smp.simplify(normalised_wf2(t) * smp.conjugate(normalised_wf1(t))),
                                     (t, self.lowerLimit, self.upperLimit), meijerg=True)
            return helper
        wavefunction = lambda x: other * self.wavefunction(x)
        return self.update_state(wavefunction)
    
    def __add__(self, other):
        if isinstance(other, QuantumState):
            def new_wavefunction(t):
                return self.wavefunction(t) + other.wavefunction(t)
            return self.update_state(new_wavefunction)
        raise TypeError('Can only add QuantumStates together')
    
    def update_state(self, wavefunction : smp.Function):
        return QuantumState(wavefunction, self.lowerLimit, self.upperLimit)
    
    def expectation_value(self, operator):
        state = QuantumState(self.normalised_wavefunction(), self.lowerLimit, self.upperLimit)
        def helper(t : smp.Symbol):
            return smp.integrate(smp.simplify(smp.conjugate(state.wavefunction(t)) * operator.apply(state).wavefunction(t)),
                                 (t, state.lowerLimit, state.upperLimit), meijerg=True)
        return helper

    #Private Methods
    def __normalise_wavefunction_exact(wavefunction : smp.Function, lowerLimit, upperLimit) -> smp.Function:
        def helper(t : smp.Symbol):
            value = smp.integrate(smp.simplify(wavefunction(t) * smp.conjugate(wavefunction(t))),
                                  (t, lowerLimit, upperLimit), meijerg=True)
            value_mod_squared = value * smp.conjugate(value)
            if value_mod_squared == smp.oo or value_mod_squared == 0:
                raise ValueError('Wavefunction is non-normalisable')
            normalised_wavefunction = 1/smp.sqrt(value) * wavefunction(t)
            return normalised_wavefunction
        return helper
    
    def __normalise_wavefunction_numer(wavefunction : smp.Function, lowerLimit, upperLimit) -> smp.Function:
        def helper(t : smp.Symbol):
            value = smp.Integral(wavefunction(t) * smp.conjugate(wavefunction(t)),
                                  (t, lowerLimit, upperLimit)).evalf()
            value_mod_squared = value * smp.conjugate(value)
            if value_mod_squared == smp.oo or value_mod_squared == 0:
                raise ValueError('Wavefunction is non-normalisable')
            normalised_wavefunction = 1/smp.sqrt(value) * wavefunction(t)
            return normalised_wavefunction
        return helper