import tomli_w
import os
import itertools
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Creates config files for simulations')

parser.add_argument('--model', type=str)                        # Which Hamiltonian to use
parser.add_argument('--length', type=int, nargs='+')            # System sizes
parser.add_argument('--disorder', type=float, nargs='+')        # Disorder values
parser.add_argument('--reals', type=int, nargs='?')             # Number of realizations
parser.add_argument('--J', type=float)                          # Coupling/Hopping parameter
parser.add_argument('--G1L', type=float)                        # Bath coupling to sigma_1^-
parser.add_argument('--G2L', type=float)                        # Bath coupling to sigma_1^+
parser.add_argument('--G1R', type=float, nargs='?', default=0)  # Bath coupling to sigma_L^-
parser.add_argument('--G2R', type=float, nargs='?', default=0)  # Bath coupling to sigma_L^+
parser.add_argument('--vectors', action='store_true')
args = parser.parse_args()

if not os.path.exists("configs"):
    os.makedirs("configs")

if args.model == "anderson":
    count = 1
    for x in itertools.product(args.length, args.disorder):
        for i in range(args.reals):

            seed = np.random.randint(1000, 100000000000)
            config = {"model": args.model, "length": x[0], "disorder": x[1], "J": args.J,
                      "gamma_1L": args.G1L, "gamma_2L": args.G2L, "gamma_1R": args.G1R, "gamma_2R": args.G2R,
                      "vectors": args.vectors, "seed": seed}
            with open(f"./configs/conf{count}.toml", "wb") as f:
                tomli_w.dump(config, f)
            count += 1

elif args.model == "aubry_andre":
    count = 1
    for x in itertools.product(args.length, args.disorder):
        config = {"model": args.model, "length": x[0], "disorder": x[1], "J": args.J,
                  "gamma_1L": args.G1L, "gamma_2L": args.G2L, "gamma_1R": args.G1R, "gamma_2R": args.G2R,
                  "vectors": args.vectors}
        with open(f"./configs/conf{count}.toml", "wb") as f:
            tomli_w.dump(config, f)
        count += 1

for L in args.length:
    for W in args.disorder:
        if not os.path.exists(f"data/L{L}/W{W}"):
            os.makedirs(f"data/L{L}/W{W}")
