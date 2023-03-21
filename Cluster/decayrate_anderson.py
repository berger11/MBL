import numpy as np
import mpmath
from flint import acb_mat, ctx
from precision_functions import *

ctx.dps = 200
def compute_gap(L, W, seed, J, gamma_1L, gamma_2L, gamma_1R=0, gamma_2R=0):

    filename = f'L{L}_W{W}_seed{seed}_J{J}_G1L{gamma_1L}_G2L{gamma_2L}_G1R{gamma_1R}_G2R{gamma_2R}'

    J = [J]
    M = build_M_matrix_edge_baths_prec(L, gamma_1L, gamma_2L, gamma_1R, gamma_2R)
    h = [-W/2, W/2]  # Convention of a factor 1/2 in front of field term
    G = build_majorana_hamiltonian_disorderedXY_prec(L, J, J, h)
    A = build_shape_matrix_prec(G, M)

    B = acb_mat(A)

    count = 0
    while count < 10:
        try:
            w = B.eig()
            break
        except ValueError:
            ctx.dps += 10
            print('Precision up')
            count += 1

    if count == 10:
        raise Exception('Matrix could not be diagonalized')

    gap = np.abs(np.real(w[0]))
    for j in range(len(w)):
        if 0 < np.real(w[j]) < gap:
            gap = np.real(w[j])
    gap = 2 * gap
    log_gap = float(gap.log())
    value = np.array([log_gap, L, W, seed, J, gamma_1L, gamma_2L, gamma_1R, gamma_2R])

    np.save(filename, value)
