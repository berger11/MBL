import numpy as np
from scipy.stats import ortho_group

def build_hamiltonian_anderson(L, t, W):
    # Returns the Hamiltonian, an LxL numpy array, of the one-dimensional Anderson model with L sites.
    # t is the hopping amplitude, W is the disorder strength

    H = np.zeros((L, L))

    for i in range(L):
        for j in range(L):

            if i == j+1 or i+1 == j:
                H[i, j] = t
            elif i == j:
                h = np.random.rand()
                h = (h - 0.5) * 2*W
                H[i, j] = h

    return H


def build_hamiltonian_aubry_andre(L, t, lambd, phase=0):
    # Returns the Hamiltonian, an LxL numpy array, of the one-dimensional Aubry-André model with L sites.
    # t is the hopping amplitude, lambd is the disorder strength

    H = np.zeros((L, L))
    x = 0.5*(np.sqrt(5) - 1)

    for i in range(L):
        for j in range(L):

            if i == j + 1 or i + 1 == j:
                H[i, j] = t
            elif i == j:
                H[i, j] = lambd * np.cos(2 * np.pi * x * (i+1) + phase)

    return H


def convert_to_majorana(H):
    # Changes the basis of the quadratic Hamiltonian H from canonical fermions to Majorana fermions
    # Accepts an LxL numpy array H and returns a 2Lx2L numpy array G.

    L = len(H)

    G = np.zeros((2*L, 2*L), dtype=np.complex_)

    for k in range(1, 2*L + 1):
        for l in range(1, 2*L + 1):

            if k%2 != 0 and l%2 != 0:

                G[k-1, l-1] = 0.25 * H[int(0.5*(k+1)) - 1, int(0.5*(l+1)) - 1]

            elif k%2 != 0 and l%2 == 0:

                G[k-1, l-1] = -0.25j * H[int(0.5*(k+1)) - 1, int(0.5*l) - 1]

            elif k%2 == 0 and l%2 != 0:

                G[k-1, l-1] = 0.25j * H[int(0.5*k) - 1, int(0.5*(l+1)) - 1]

            elif k%2 == 0 and l%2 == 0:

                G[k-1, l-1] = 0.25 * H[int(0.5*k) - 1, int(0.5*l) - 1]

    G = 0.5*(G - G.T) #Antisymmetrizes
    return G


def build_shape_matrix(G, M):
    # Builds the shape matrix as defined in eq. (27) in Prosen
    # Takes as argument a 2L x 2L Hamiltonian G in the Majorana basis and the 2L x 2L M matrix as defined in eq. (23).

    L = int(len(G) / 2)

    A = np.zeros((4*L, 4*L), dtype=np.complex_)

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


def generate_binary_strings(m):
    # Generates all binary sequences of length m, and returns them as a list of strings

    bi = bin(2**m - 1)

    line = [bin(i)[2:] for i in range(int(str(bi), 2) + 1)]

    maxlength = len(line[-1])

    for i in range(len(line)):
        if len(line[i]) < maxlength:
            for j in range(maxlength - len(line[i])):
                line[i] = "0" + line[i]

    return line


def compute_spectrum(beta):
    # Computes the spectrum of the Lindbladian from the rapidities of the shape matrix, using eq (42) in Prosen
    # Takes as argument a numpy array of 2*L rapidities
    # Returns a list of 4**L Lindblad eigenvalues

    spectrum = np.array([])

    binary = generate_binary_strings(2**(len(beta))-1)

    for bi in binary:
        bin_list = np.array(list(bi))
        prod = bin_list * beta
        eigval = -2 * np.sum(prod)
        np.append(spectrum, eigval)

    return spectrum


def build_majorana_hamiltonian_homogenousXY(L, J_x, J_y, h):
    # Builds the Hamiltonian in the Majorana basis of the homogenous XY model (Prosen sec. 4.1)
    # Returns a 2L x 2L numpy array

    G = np.zeros((2*L, 2*L), dtype=np.complex_)

    for m in range(1, L+1):

        G[2*m-2, 2*m-1] -= 1j * h
        if m != L:
            G[2*m-1, 2*m] -= 1j * J_x
            G[2*m-2, 2*m + 1] += 1j * J_y

    G = 0.5 * (G - G.T)

    return G


def build_majorana_hamiltonian_disorderedXY(L, J_x, J_y, h, seed=None):
    # Builds the Hamiltonian in the Majorana basis of the disordered XY model (Prosen eq. (52))
    # J_x, J_y and h are numpy arrays with either one or two entries
    # If the array only has one entry, then the parameter will be homogenous
    # If the array has two entries, the parameter will be drawn from a uniform distribution with lower and upper bound given by the first and second entry respectively
    # Returns a 2L x 2L numpy array

    G = np.zeros((2 * L, 2 * L), dtype=np.complex_)

    if seed is not None:
        np.random.seed(seed)

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


def build_M_matrix_edge_baths(L, gamma_1L, gamma_2L, gamma_1R, gamma_2R):
    # Builds the M matrix as defined in eq. 23 in Prosen, where baths are coupled to the two edges of a chain
    # L is the system size, gamma are the bath coupling parameters as defined in eq. 55
    # Returns a 2L x 2L numpy array

    M = np.zeros((2*L, 2*L), dtype=np.complex_)

    M[0, 0] = 0.25 * gamma_1L + 0.25 * gamma_2L
    M[0, 1] = 0.25j * gamma_1L - 0.25j * gamma_2L
    M[1, 0] = -0.25j * gamma_1L + 0.25j * gamma_2L
    M[1, 1] = 0.25 * gamma_1L + 0.25 * gamma_2L

    M[-2, -2] = 0.25 * gamma_1R + 0.25 * gamma_2R
    M[-2, -1] = 0.25j * gamma_1R - 0.25j * gamma_2R
    M[-1, -2] = -0.25j * gamma_1R + 0.25j * gamma_2R
    M[-1, -1] = 0.25 * gamma_1R + 0.25 * gamma_2R

    return M


def build_M_matrix_disordered_bath_left(L, trD=2):
    # Builds the M matrix as defined eq. 23 in Prosen, but using the disordered coupling as in Huse.
    # L is the system size. The trace of the D matrix in Huse can be chosen, default is 2
    # Instead of (Sigma^X, Sigma^Y, Sigma^Z), the bath operators are chosen as (Sigma^+, Sigma^-).
    # Returns a 2L x 2L numpy array

    M = np.zeros((2 * L, 2 * L), dtype=np.complex_)

    D_1 = np.random.uniform(0, 1)
    D_2 = np.random.uniform(0, 1)

    D_1 = trD * D_1 / (D_1 + D_2)
    D_2 = trD * D_2 / (D_1 + D_2)

    U = ortho_group.rvs(2)

    l_1_1 = np.sqrt(D_1 / 2) * (U[0, 0] + U[0, 1])
    l_1_2 = 1j * np.sqrt(D_1 / 2) * (U[0, 1] - U[0, 0])
    l_2_1 = np.sqrt(D_2 / 2) * (U[1, 0] + U[1, 1])
    l_2_2 = 1j * np.sqrt(D_2/2) * (U[1, 1] - U[1, 0])

    M[0, 0] = l_1_1 * np.conjugate(l_1_1) + l_2_1 * np.conjugate(l_2_1)
    M[1, 0] = l_1_2 * np.conjugate(l_1_1) + l_2_2 * np.conjugate(l_2_1)
    M[0, 1] = l_1_1 * np.conjugate(l_1_2) + l_2_1 * np.conjugate(l_2_2)
    M[1, 1] = l_1_2 * np.conjugate(l_1_2) + l_2_2 * np.conjugate(l_2_2)

    return M












