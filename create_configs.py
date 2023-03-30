import tomli_w
import os
import itertools
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Creates config files for simulations')

parser.add_argument('--length', type=int, nargs='+')
parser.add_argument('--disorder', type=float, nargs='+')
parser.add_argument('--reals', type=int)
parser.add_argument('--model', type=str)
args = parser.parse_args()

L_list = args.length
W_list = args.disorder
reals = args.reals

if not os.path.exists("configs"):
    os.makedirs("configs")
if args.model == "anderson":
    count = 1
    for x in itertools.product(L_list, W_list):
        for i in range(reals):

            seed = np.random.randint(1000, 100000000000)
            config = {"length": x[0], "disorder": x[1], "seed": seed, "J": 0.25,
                      "gamma_1L": 0.5, "gamma_2L": 0.5, "gamma_1R": 0, "gamma_2R": 0}
            with open(f"./configs/conf{count}.toml", "wb") as f:
                tomli_w.dump(config, f)
            count += 1
elif args.model == "aubry_andre":
    for x in itertools.product(L_list, W_list):
        config = {"length": x[0], "disorder": x[1], "J": 0.25,
                  "gamma_1L": 0.5, "gamma_2L": 0.5, "gamma_1R": 0, "gamma_2R": 0}
for L in L_list:
    for W in W_list:
        if not os.path.exists(f"data/L{L}/W{W}"):
            os.makedirs(f"data/L{L}/W{W}")
