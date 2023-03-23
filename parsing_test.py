import argparse
from functions.decayrate_anderson import compute_gap
import numpy as np
parser = argparse.ArgumentParser(description='Does stuff')

parser.add_argument('--seed', type=int)

args = parser.parse_args()

filename, value = compute_gap(5, 1, args.seed, 0.25, 0.5, 0.5)

np.save(filename, value)