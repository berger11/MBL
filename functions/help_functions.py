import numpy as np
#from scipy.stats import unitary_group
from scipy.stats import ortho_group
from scipy import sparse

def Build_Pauli_Mats(L, k):
    # Returns the matrix representations of the three Pauli matrices at site k in a chain of size dim
    # leftmost site has index 0

    X = sparse.csr_matrix(np.array([[0, 1], [1, 0]]))
    Y = sparse.csr_matrix(np.array([[0, 0 - 1j], [0 + 1j, 0]]))
    Z = sparse.csr_matrix(np.array([[1, 0], [0, -1]]))
    Id = sparse.csr_matrix(np.array([[1, 0], [0, 1]]))

    if k == 0:

        X_full = X
        Y_full = Y
        Z_full = Z

    else:

        X_full = sparse.csr_matrix(np.array([[1, 0], [0, 1]]))
        Y_full = sparse.csr_matrix(np.array([[1, 0], [0, 1]]))
        Z_full = sparse.csr_matrix(np.array([[1, 0], [0, 1]]))

    for i in range(1, L):

        if i == k:

            X_full = sparse.kron(X_full, X)
            Y_full = sparse.kron(Y_full, Y)
            Z_full = sparse.kron(Z_full, Z)

        else:

            X_full = sparse.kron(X_full, Id)
            Y_full = sparse.kron(Y_full, Id)
            Z_full = sparse.kron(Z_full, Id)

    return X_full, Y_full, Z_full


def build_Kossakowski_matrix():
    # Returns the Kossakowski matrix of size 3 x 3 as a numpy array

    eig1 = np.random.uniform(low=0, high=1)
    eig2 = np.random.uniform(low=0, high=1)
    eig3 = np.random.uniform(low=0, high=1)

    D = np.array([[eig1, 0, 0],
                  [0, eig2, 0],
                  [0, 0, eig3]])

    trD = np.trace(D)

    D = 0.5 * D/trD

    U = ortho_group.rvs(3)

    K = U.T @ D @ U

    return K

def Build_Hamiltonian_XXX(L, W, X_list, Y_list, Z_list):
    # Builds XXX Hamiltonian of a chain of size L. Returns 2^L x 2^L matrix representation of Hamiltonian.
    # W is the disorder strength
    # Sigma_list is the list of Pauli matrices for L sites
    # returns the hamiltonian, as well as the interaction term for testing purposes

    interaction_term = np.zeros((2**L, 2**L))

    for i in range(L-1):
        X_i = X_list[i]
        Y_i = Y_list[i]
        Z_i = Z_list[i]

        X_i_plus_1 = X_list[i + 1]
        Y_i_plus_1 = Y_list[i + 1]
        Z_i_plus_1 = Z_list[i + 1]
        #print(Z_i @ Z_i_plus_1)
        ith_term = 0.25 * (X_i @ X_i_plus_1 + Y_i @ Y_i_plus_1 + Z_i @ Z_i_plus_1)
        interaction_term = interaction_term + ith_term

    field_term = np.zeros((2**L, 2**L))

    for j in range(L):
        field = np.random.uniform(low=-W, high=W)
        Z_j = Z_list[j]
        field_term = field_term + 0.5 * field * Z_j

    hamiltonian = interaction_term + field_term

    return hamiltonian, interaction_term

def build_lindblad_matrix(L, H, K, X_0, Y_0, Z_0):
    # Returns the matrix representation of the Lindblad superoperator of size 4**L x 4**L as a numpy array
    # L is the number of sites, H is the 2**L x 2**L Hamiltonian, K is the 3x3 Kossakowski matrix.
    # X_0, Y_0 and Z_0 are the 2**L x 2**L Pauli matrices at the leftmost site

    Id = np.identity(2 ** L)

    #X_0 = X_list[0]
    #Y_0 = Y_list[0]
    #Z_0 = Z_list[0]

    commutator_term = -1j * np.kron(H, Id) + 1j * np.kron(Id, H)

    term_00 = K[0, 0] * (sparse.kron(X_0, X_0) - 0.5 * sparse.kron(X_0 @ X_0, Id) - 0.5 * sparse.kron(Id, (X_0 @ X_0).T))
    term_01 = K[0, 1] * (sparse.kron(X_0, Y_0.T) - 0.5 * sparse.kron(Y_0 @ X_0, Id) - 0.5 * sparse.kron(Id, (Y_0 @ X_0).T))
    term_02 = K[0, 2] * (sparse.kron(X_0, Z_0) - 0.5 * sparse.kron(Z_0 @ X_0, Id) - 0.5 * sparse.kron(Id, (Z_0 @ X_0).T))
    term_10 = K[1, 0] * (sparse.kron(Y_0, X_0) - 0.5 * sparse.kron(X_0 @ Y_0, Id) - 0.5 * sparse.kron(Id, (X_0 @ Y_0).T))
    term_11 = K[1, 1] * (sparse.kron(Y_0, Y_0.T) - 0.5 * sparse.kron(Y_0 @ Y_0, Id) - 0.5 * sparse.kron(Id, (Y_0 @ Y_0).T))
    term_12 = K[1, 2] * (sparse.kron(Y_0, Z_0) - 0.5 * sparse.kron(Z_0 @ Y_0, Id) - 0.5 * sparse.kron(Id, (Z_0 @ Y_0).T))
    term_20 = K[2, 0] * (sparse.kron(Z_0, X_0) - 0.5 * sparse.kron(X_0 @ Z_0, Id) - 0.5 * sparse.kron(Id, (X_0 @ Z_0).T))
    term_21 = K[2, 1] * (sparse.kron(Z_0, Y_0.T) - 0.5 * sparse.kron(Y_0 @ Z_0, Id) - 0.5 * sparse.kron(Id, (Y_0 @ Z_0).T))
    term_22 = K[2, 2] * (sparse.kron(Z_0, Z_0) - 0.5 * sparse.kron(Z_0 @ Z_0, Id) - 0.5 * sparse.kron(Id, (Z_0 @ Z_0).T))

    sum_term = term_00 + term_01 + term_02 + term_10 + term_11 + term_12 + term_20 + term_21 + term_22

    lindblad = commutator_term + sum_term

    return lindblad



