import tomli_w
import os
import itertools
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Creates config files for simulations')

parser.add_argument('--sim_name', type=str)                       # Name of directory where data is saved
parser.add_argument('--model', type=str)                          # Which Hamiltonian to use
parser.add_argument('--length', type=int, nargs='+')              # System sizes
parser.add_argument('--disorder', type=float, nargs='+')          # Disorder values
parser.add_argument('--reals', type=int, nargs='+')               # Number of realizations
parser.add_argument('--t', type=float)                            # Hopping parameter
parser.add_argument('--G1L', type=float, nargs='?', default=0.5)  # Bath coupling to sigma_1^-
parser.add_argument('--G2L', type=float, nargs='?', default=0.5)  # Bath coupling to sigma_1^+
parser.add_argument('--G1R', type=float, nargs='?', default=0)    # Bath coupling to sigma_L^-
parser.add_argument('--G2R', type=float, nargs='?', default=0)    # Bath coupling to sigma_L^+
parser.add_argument('--vectors', action='store_true')             # Whether to compute eigenvectors or not
args = parser.parse_args()

if not os.path.exists("configs"):
    os.makedirs("configs")


for L_ind, L in enumerate(args.length):
    for W in args.disorder:
        config = {}
        for arg in vars(args):
            config[arg] = getattr(args, arg)

for L in args.length:
    for W in args.disorder:
        if not os.path.exists(f"{args.sim_name}/L{L}/W{W}"):
            os.makedirs(f"{args.sim_name}/L{L}/W{W}")
