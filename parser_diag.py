import argparse
from functions.model_diagonalization import eigen_anderson, eigen_aubryandre
import numpy as np
import tomli

parser = argparse.ArgumentParser(description='Does stuff')

parser.add_argument('--index', type=int)
args = parser.parse_args()

with open(f"./configs/conf{args.index}.toml", "rb") as f:
    conf = tomli.load(f)

if conf["model"] == "anderson":

    L = conf["length"]
    W = conf["disorder"]
    seed = conf["seed"]
    J = conf["J"]
    gamma_1L = conf["gamma_1L"]
    gamma_2L = conf["gamma_2L"]
    gamma_1R = conf["gamma_1R"]
    gamma_2R = conf["gamma_2R"]
    eigvecs = conf["vectors"]

    value = eigen_anderson(L, W, seed, J, gamma_1L, gamma_2L, gamma_1R, gamma_2R, eigvecs)

    np.save(f"data_anderson/L{L}/W{W}/anderson_ind{args.index}_seed{seed}", value)
    print(f"Saved data anderson L{L} W{W}")

elif conf["model"] == "aubry_andre":

    L = conf["length"]
    W = conf["disorder"]
    t = conf["J"]
    gamma_1L = conf["gamma_1L"]
    gamma_2L = conf["gamma_2L"]
    gamma_1R = conf["gamma_1R"]
    gamma_2R = conf["gamma_2R"]
    eigvecs = conf["vectors"]

    value = eigen_aubryandre(L, W, t, gamma_1L, gamma_2L, gamma_1R, gamma_2R, eigvecs)

    np.save(f"data_aubry_andre/L{L}/W{W}/aubryandre_ind{args.index}", value)
    print(f"Saved data aubry andre L{L} W{W}")
