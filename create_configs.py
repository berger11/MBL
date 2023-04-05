import tomli_w
import os
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

if not os.path.exists(f'configs/{args.sim_name}_configs'):
    os.makedirs(f'configs/{args.sim_name}_configs')

count = 0
for L_ind, L in enumerate(args.length):
    for W in args.disorder:
        config = {'length': L, 'disorder': W, 'reals': args.reals[L_ind]}
        for arg in vars(args):
            if arg != 'length' and arg != 'disorder' and arg != 'reals':
                config[arg] = getattr(args, arg)
        with open(f'configs/{args.sim_name}_configs/{args.sim_name}_conf{count}.toml', 'wb') as f:
            tomli_w.dump(config, f)
        count += 1
            
for L in args.length:
    for W in args.disorder:
        if not os.path.exists(f'data/{args.sim_name}_data/L{L}/W{W}'):
            os.makedirs(f'data/{args.sim_name}_data/L{L}/W{W}')
