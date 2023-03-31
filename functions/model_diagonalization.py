from flint import acb_mat, ctx
import numpy as np
from functions.precision_functions import *
from functions.tq_help_functions import build_hamiltonian_aubry_andre, convert_to_majorana
from functions.flint_help_functions import _spectrum_to_string
import mpmath

ctx.dps = 200


def eigen_anderson(L, W, seed, J, gamma_1L, gamma_2L, gamma_1R, gamma_2R, comp_vecs):

    #filename = f'./data/L{L}_W{W}_seed{seed}_J{J}_G1L{gamma_1L}_G2L{gamma_2L}_G1R{gamma_1R}_G2R{gamma_2R}'

    J = [J]
    M = build_M_matrix_edge_baths_prec(L, gamma_1L, gamma_2L, gamma_1R, gamma_2R)
    h = [-W/2, W/2]  # Convention of a factor 1/2 in front of field term
    G = build_majorana_hamiltonian_disorderedXY_prec(L, J, J, h)
    A = build_shape_matrix_prec(G, M)

    B = acb_mat(A)

    count = 0
    while count < 10:
        try:
            eigens = B.eig(right=comp_vecs)
            break
        except ValueError:
            ctx.dps += 10
            print('Precision up')
            count += 1

    if count == 10:
        raise Exception('Matrix could not be diagonalized')

    """gap = np.abs(np.real(w[0]))
    for j in range(len(w)):
        if 0 < np.real(w[j]) < gap:
            gap = np.real(w[j])
    gap = 2 * gap
    log_gap = float(gap.log())"""

    eigens = _spectrum_to_string(L, eigens, comp_vecs)

    value = np.array([eigens, L, W, seed, J[0], gamma_1L, gamma_2L, gamma_1R, gamma_2R, comp_vecs], dtype=object)

    return value


def eigen_aubryandre(L, lambd, t, gamma_1L, gamma_2L, gamma_1R, gamma_2R, comp_vecs):

    #filename = f'./data/aub_L{L}_lambd{lambd}_t{t}_G1L{gamma_1L}_G2L{gamma_2L}_G1R{gamma_1R}_G2R{gamma_2R}'

    M = build_M_matrix_edge_baths_prec(L, gamma_1L, gamma_2L, gamma_1R, gamma_2R)
    H = build_hamiltonian_aubry_andre(L, t, lambd)
    G = convert_to_majorana(H)
    G = mpmath.matrix(G)
    A = build_shape_matrix_prec(G, M)

    B = acb_mat(A)

    count = 0
    while count < 10:
        try:
            eigens = B.eig(right=comp_vecs)
            break
        except ValueError:
            ctx.dps += 10
            print('Precision up')
            count += 1

    if count == 10:
        raise Exception('Matrix could not be diagonalized')

    """gap = np.abs(np.real(w[0]))
    for j in range(len(w)):
        if 0 < np.real(w[j]) < gap:
            gap = np.real(w[j])
    gap = 2 * gap
    log_gap = float(gap.log())"""

    eigens = _spectrum_to_string(L, eigens, comp_vecs)

    value = np.array([eigens, L, lambd, t, gamma_1L, gamma_2L, gamma_1R, gamma_2R, comp_vecs], dtype=object)

    return value