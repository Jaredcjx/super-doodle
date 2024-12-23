from sympy.physics.quantum.constants import hbar
import sympy as smp
import numpy as np
import matplotlib.pyplot as plt
import QuantumBase as qb

Hamiltonian : qb.Operator
State : qb.PositionQuantumState
m, w, x = smp.symbols('m w x', real = True, positive = True)
P, X = qb.MomentumOperator(), qb.PositionOperator()
'''
Creating the Hamiltonian for Simple Harmonic Motion: H = hw(N + 1/2 I), N = a_dagger * a
Hamiltonian_Momentum = P ** 2 * (1/(2*m))
Hamiltonian_Potential = X ** 2 * (m * w ** 2 / 2)
'''
a_dagger = X * smp.sqrt(m * w / (2 * hbar)) - P * (smp.I / (smp.sqrt(2 * m * w * hbar)))
a = X * smp.sqrt(m * w / (2 * hbar)) + P * (smp.I / (smp.sqrt(2 * m * w * hbar)))
N = a_dagger * a
Hamiltonian = (N + qb.IdentityOperator() * (1/2)) * hbar * w

def get_eigenvalue_En(n : int):
    return hbar * w * (n + 1/2)

def create_eigenstate_En(n: int, state_base = None):
    if state_base == None:
        solutions_0 = solve_base_state()
        psi_0 = get_solution(solutions_0)
        state_0 = qb.PositionQuantumState(psi_0)
        state_n = prepare_state_n(state_0, n)
    else:
        wavefunction = state_base.wavefunction(x)
        state_base = state_base.update_state(lambda y: wavefunction.subs(x, y))
        state_n = prepare_state_n(state_base, n)
    return state_n

def solve_base_state():
    # Solve H|n> = E_n|n>
    def helper(t : smp.Symbol):
        psi_0 = smp.Function('psi_0')
        state_0 = qb.PositionQuantumState(psi_0)
        solutions_0 = smp.dsolve(smp.Eq(a.apply(state_0).wavefunction(t), 0), psi_0(t))
        return solutions_0
    return helper

def get_solution(solutions_0):
    def helper(t : smp.Symbol):
        solution = solutions_0(t)
        res = solution.rhs.xreplace({s: 1 for s in solution.rhs.atoms(smp.Symbol) if 'C' in s.name})
        return res
    return helper

def prepare_state_n(state_base, n):
    a_dagger_n = a_dagger ** n
    state_n = a_dagger_n.apply(state_base)
    return state_n

def time_evolution(state, time, num_terms = None):
    eigenstate_eigenvalue_list = [(create_eigenstate_En(0),
                                   smp.exp(-smp.I * get_eigenvalue_En(0) * time / hbar))]
    if num_terms == None:
        num_terms = 10
    for i in range(1, num_terms):
        eigenstate_eigenvalue_list.append((create_eigenstate_En(1, eigenstate_eigenvalue_list[-1][0]),
                                           smp.exp(-smp.I * get_eigenvalue_En(i) * time / hbar)))
    final_State = None
    for eigenstate, eigenvalue in eigenstate_eigenvalue_list:
        c_n = state * eigenstate
        c_n = c_n(x)
        if final_State == None:
            final_State = eigenstate * (c_n * eigenvalue)
        else:
            final_State += eigenstate * (c_n * eigenvalue)
    return final_State

def graph_plot(initial_wavefunction, final_wavefunction, normalise_numerical : bool = None):
    x_vals = np.linspace(-5, 5, 1000)
    f = lambda x: final_wavefunction(x).subs(hbar, 1).subs(m, 1).subs(w, 1)
    g = lambda x: initial_wavefunction(x).subs(hbar, 1).subs(m, 1).subs(w, 1)
    
    if normalise_numerical:
        initial_state = qb.QuantumState(g)
        final_state = qb.QuantumState(f)
        initial_wavefunction_func = smp.lambdify(x, initial_state.normalised_wavefunction(numerical=True)(x), 'numpy')
        final_wavefunction_func = smp.lambdify(x, final_state.normalised_wavefunction(numerical=True)(x), 'numpy')
    else:
        initial_wavefunction_func = smp.lambdify(x, g(x), 'numpy')
        final_wavefunction_func = smp.lambdify(x, f(x), 'numpy')

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, initial_wavefunction_func(x_vals), label="Initial Wavefunction: $\\psi(x)$")
    plt.plot(x_vals, final_wavefunction_func(x_vals), label="Final Wavefunction: $\\psi_{final}(x)$", linestyle='--')
    plt.title("Initial and Final Wavefunctions")
    plt.xlabel("x")
    plt.ylabel("$\\psi(x)$")
    plt.legend()
    plt.grid()
    plt.show()

def main():
    '''
    Creating a sample wavefunction: psi(x) = exp(-x**2)
    '''
    Wavefunction = lambda x: smp.exp(-x**2)
    initial_state = qb.PositionQuantumState(Wavefunction)
    
    final_state = time_evolution(state=initial_state, time=2 * smp.pi, num_terms=10)

    initial_state_momentum = qb.MomentumQuantumState(initial_state.get_momentum_basis_wavefunction())
    final_state_momentum = qb.MomentumQuantumState(final_state.get_momentum_basis_wavefunction())
    
    graph_plot(initial_state.wavefunction, final_state.wavefunction, normalise_numerical=True)
    graph_plot(initial_state_momentum.wavefunction, final_state_momentum.wavefunction, normalise_numerical=True)

if __name__ == "__main__":
    main()