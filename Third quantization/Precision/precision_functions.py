import numpy as np
import mpmath

# In this program I recreate some of the functions of tq_help_functions.py, but with mpmath

def build_majorana_hamiltonian_disorderedXY_prec(L, J_x, J_y, h):
    # Builds the Hamiltonian in the Majorana basis of the disordered XY model (Prosen eq. (52))
    # J_x, J_y and h are lists with either one or two entries
    # If the array only has one entry, then the parameter will be homogenous
    # If the array has two entries, the parameter will be drawn from a uniform distribution with lower and upper bound given by the first and second entry respectively
    # Returns a 2L x 2L mpmath matrix

    #G = np.zeros((2 * L, 2 * L), dtype=np.complex_)
    G = mpmath.matrix(2*L)
    for m in range(1, L+1):
        if len(h) == 1:
            G[2*m-2, 2*m-1] -= 1j * h[0]
        elif len(h) == 2:
            G[2*m-2, 2*m-1] -= 1j * np.random.uniform(low=h[0], high=h[1])

        if m != L:
            if len(J_x) == 1:
                G[2*m-1, 2*m] -= 1j * J_x[0]
            elif len(J_x) == 2:
                G[2*m-1, 2*m] -= 1j * np.random.uniform(low=J_x[0], high=J_x[1])

            if len(J_y) == 1:
                G[2*m-2, 2*m + 1] += 1j * J_y[0]
            elif len(J_y) == 2:
                G[2*m-2, 2*m+1] += 1j * np.random.uniform(low=J_y[0], high=J_y[1])

    G = 0.5 * (G - G.T)

    return G


def build_M_matrix_edge_baths_prec(L, gamma_1L, gamma_2L, gamma_1R, gamma_2R):
    # Builds the M matrix as defined in eq. 23 in Prosen, where baths are coupled to the two edges of a chain
    # L is the system size, gamma are the bath coupling parameters as defined in eq. 55
    # Returns a 2L x 2L mpmath matrix

    #M = np.zeros((2*L, 2*L), dtype=np.complex_)
    M = mpmath.matrix(2*L)

    M[0, 0] = 0.25 * gamma_1L + 0.25 * gamma_2L
    M[0, 1] = 0.25j * gamma_1L - 0.25j * gamma_2L
    M[1, 0] = -0.25j * gamma_1L + 0.25j * gamma_2L
    M[1, 1] = 0.25 * gamma_1L + 0.25 * gamma_2L

    M[-2, -2] = 0.25 * gamma_1R + 0.25 * gamma_2R
    M[-2, -1] = 0.25j * gamma_1R - 0.25j * gamma_2R
    M[-1, -2] = -0.25j * gamma_1R + 0.25j * gamma_2R
    M[-1, -1] = 0.25 * gamma_1R + 0.25 * gamma_2R

    return M

def build_shape_matrix_prec(G, M):
    # Builds the shape matrix as defined in eq. (27) in Prosen
    # Takes as argument a 2L x 2L Hamiltonian G in the Majorana basis and the 2L x 2L M matrix as defined in eq. (23).

    L = int(len(G) / 2)

    #A = np.zeros((4*L, 4*L), dtype=np.complex_)
    A = mpmath.matrix(4*L)

    for n in range(1, 4*L + 1):
        for m in range(1, 4*L + 1):

            if n%2 != 0 and m%2 != 0:

                A[n-1, m-1] = -2j * G[int(0.5*(n+1)) - 1, int(0.5*(m+1)) - 1] - M[int(0.5*(n+1)) - 1, int(0.5*(m+1)) - 1] + M[int(0.5*(m+1)) - 1, int(0.5*(n+1)) - 1]

            elif n%2 != 0 and m%2 == 0:

                A[n-1, m-1] = 2j * M[int(0.5*m) - 1, int(0.5*(n+1)) - 1]

            elif n%2 == 0 and m%2 != 0:

                A[n-1, m-1] = -2j * M[int(0.5*n) - 1, int(0.5*(m+1)) - 1]

            elif n%2 == 0 and m%2 == 0:

                A[n-1, m-1] = -2j * G[int(0.5*n) - 1, int(0.5*m) - 1] + M[int(0.5*n) - 1, int(0.5*m) - 1] - M[int(0.5*m) - 1, int(0.5*n) - 1]

    return A