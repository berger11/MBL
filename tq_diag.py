import argparse
import tomli
from functions.tq_help_functions import (build_hamiltonian_aubry_andre, build_hamiltonian_anderson, convert_to_majorana,
                               build_M_matrix_edge_baths, build_shape_matrix)
from flint import ctx, acb_mat, arb
from mpmath import matrix as mpmath_matrix
from functions.flint_help_functions import _spectrum_to_string
import numpy as np

parser = argparse.ArgumentParser(description='Diagonalizes shape matrix of a single Hamiltonian realization')

parser.add_argument('--sim_name', type=str)  # Name of simulation
parser.add_argument('--task_id', type=int)   # Task ID as given from slurm job array
parser.add_argument('--config_index', type=int)  # Index of the configuration file
args = parser.parse_args()

with open(f'configs/{args.sim_name}_configs/{args.sim_name}_conf{args.config_index}.toml', 'rb') as f:
    config = tomli.load(f)

model = config['model']
L = config['length']
W = config['disorder']
t = config['t']
G1L = config['G1L']
G2L = config['G2L']
G1R = config['G1R']
G2R = config['G2R']
vectors = config['vectors']
precision = config['precision']

ctx.dps = precision

seed = np.random.randint(1, 2**31)
config['seed'] = seed

H = eval(f'build_hamiltonian_{model}({L}, {t}, {W}, seed={seed})')
G = convert_to_majorana(H)
M = build_M_matrix_edge_baths(L, G1L, G2L, G1R, G2R)
A = build_shape_matrix(G, M)
A = mpmath_matrix(A)  # Not possible to convert straight from numpy array to acb_mat
A = acb_mat(A)

count = 0
while count < 10:
    try:
        eigens = A.eig(right=vectors)
        break
    except ValueError:
        ctx.dps += 10
        print('Precision up')
        count += 1

if count == 10:
    raise Exception('Matrix could not be diagonalized')

eigens = _spectrum_to_string(L, eigens, vectors)
value = np.array([eigens, config], dtype=object)
save_path = f'data/{args.sim_name}_data/L{L}/W{W}/{args.task_id}'
np.save(save_path, value)
