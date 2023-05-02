import numpy as np
from tq_help_functions import build_hamiltonian_anderson, build_hamiltonian_aubry_andre

def compute_xi(L, E, h, t):
    # Computes the localization length of an eigenstate with energy E in a disordered 1D chain with open boundary conditions
    # L is the system size, E is the eigenstate energy, h is a 1D array of length L containing the on-site disorders, t is the hopping parameter

    M = np.array([[1, 0],
                  [0, 1]])

    for i in range(2, L):
        M_i = np.array([[(E - h[i-1])/t, -1],
                        [1, 0]])
        M = M_i @ M

    prod = M.T @ M  # Operator norm is the square root of the largest eigenvaue of M.T @ M
    #print(prod)
    w, v = np.linalg.eig(prod)

    op_norm = np.sqrt(np.max(w))
    xi = 1 / (np.log(op_norm) / (L-2))

    return xi

def get_avg_xi(L, W, t, model='anderson'):
    # Computes the average localization length over energy eigenvalues for a disordered 1D chain
    # L is the system size, W is the disorder strength, reals is the number of realizations, t is the hopping parameter

    H = eval(f'build_hamiltonian_{model}({L}, {t}, {W})')
    h = np.diagonal(H)
    w, v = np.linalg.eig(H)

    xi_avg = 0
    for E in w:

        xi = compute_xi(L, E, h, t)
        #print(xi)
        xi_avg += xi

    xi_avg = xi_avg / L

    return xi_avg
