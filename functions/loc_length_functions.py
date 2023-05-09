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
    w, v = np.linalg.eig(prod)

    op_norm = np.sqrt(np.max(w))
    xi = 1 / (np.log(op_norm) / (L-2))

    return xi

def get_avg_xi(L, W, t, model='anderson'):
    # Computes the localization length averaged over energy eigenvalues for a disordered 1D chain
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

def col_transform(M):
    # Takes as argument a 2x2 numpy array M
    # Transforms it according to Eq. 18 in https://doi.org/10.1007/BF01578242
    # Returns the transformed matrix and the normalization factors b0, b1

    b0 = np.linalg.norm(M[:, 0])
    M[:, 0] = M[:, 0] / b0
    M[:, 1] = M[:, 1] - np.dot(M[:, 0], M[:, 1]) * M[:, 0]
    b1 = np.linalg.norm(M[:, 1])
    M[:, 1] = M[:, 1] / b1

    return M, b0, b1

def compute_loc_length(L, t, E, W, model='anderson', seed=None):
    # Takes as argument system size L, hopping parameter t, energy E, disorder strength W, disorder seed seed
    # Returns final product matrix after col_transform, the c_n coefficients, the d_n coefficients

    M = np.array([[1, 0],
                  [0, 1]])
    c0 = np.zeros(L)
    c1 = np.zeros(L)
    d0 = np.zeros(L)
    d1 = np.zeros(L)

    if seed is not None:
        np.random.seed(seed)

    if model == 'aubry_andre':
        phase = np.random.uniform(low=0, high=2*np.pi)

    c0_i = 0
    c1_i = 0
    d0_i = 0
    d1_i = 0
    x = 0.5 * (np.sqrt(5) - 1)

    for i in range(L):
        if model == 'anderson':
            h_i = np.random.uniform(low=-W, high=W)
        elif model == 'aubry_andre':
            h_i = W * np.cos(2 * np.pi * x * (i+1) + phase)

        T_i = np.array([[(E - h_i) / t, -1],
                        [1, 0]])
        M = T_i @ M
        M, b0, b1 = col_transform(M)

        c0_i = c0_i + np.log(b0)
        c1_i = c1_i + np.log(b1)
        d0_i = d0_i + (np.log(b0))**2
        d1_i = d1_i + (np.log(b1))**2
        c0[i] = c0_i
        c1[i] = c1_i
        d0[i] = d0_i
        d1[i] = d1_i

    return M, c0, c1, d0, d1
