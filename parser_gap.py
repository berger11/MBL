import argparse
from functions.decayrate_anderson import compute_gap
import numpy as np
import tomli
import os

parser = argparse.ArgumentParser(description='Does stuff')

parser.add_argument('--index', type=int)

args = parser.parse_args()

with open(f"./configs/conf{args.index}.toml", "rb") as f:
    conf = tomli.load(f)

L = conf["length"]
W = conf["disorder"]
seed = conf["seed"]
J = conf["J"]
gamma_1L = conf["gamma_1L"]
gamma_2L = conf["gamma_2L"]
gamma_1R = conf["gamma_1R"]
gamma_2R = conf["gamma_2R"]

filename, value = compute_gap(L, W, seed, J, gamma_1L, gamma_2L, gamma_1R, gamma_2R)

np.save(f"data/L{L}/W{W}/ind{args.index}_seed{seed}", value)