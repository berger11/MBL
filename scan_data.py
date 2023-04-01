import numpy as np
import os
from statistics import mean
import argparse
from functions.flint_help_functions import string_to_acb
from flint import arb, acb, acb_mat, ctx
# Scans data from a given directory and computes the 10-log of the gap, and saves the values to an array

ctx.dps = 300

parser = argparse.ArgumentParser()
parser.add_argument('--directory', type=str)
parser.add_argument('--filename', type=str)
args = parser.parse_args()
L_list = []
W_list = []

for subdir, dirs, files in os.walk(args.directory):
    print(subdir)
    subdir_split = subdir.split("/")
    last = subdir_split[-1]
    if last[0] == "L":
        number = last[1:]
        if int(number) not in L_list:
            L_list.append(int(number))
    if last[0] == "W":
        number = last[1:]
        if float(number) not in W_list:
            W_list.append(float(number))


values = np.zeros((len(L_list), len(W_list)))

split_name = args.directory.split("_")
if split_name[1] == "aubry":
    for L_ind, L in enumerate(L_list):
        for W_ind, W in enumerate(W_list):
            print(f"L{L} W{W}")
            path = f"./{args.directory}/L{L}/W{W}"
            file_count = 0
            for entry in os.scandir(path):
                data = np.load(path+"/"+entry.name, allow_pickle=True)
                comp_vecs = data[-1]
                if comp_vecs:
                    spectrum = data[0][0]
                else:
                    spectrum = data[0]
                for i in range(4*L):
                    spectrum[i] = string_to_acb(spectrum[i])
                gap = np.abs(np.real(spectrum[0]))
                for j in range(len(spectrum)):
                    if 0 < np.real(spectrum[j]) < gap:
                        gap = np.real(spectrum[j])
                gap = 2 * gap
                log_gap = float(gap.log() / arb(10).log())
                values[L_ind, W_ind] = log_gap
                file_count += 1

            try:
                assert file_count == 1
            except AssertionError:
                print(f"WARNING: MORE/LESS THAN ONE FILE FOUND FOR L{L} W{W} FOR AUBRY-ANDRE MODEL")

np.save(f"npy_data/{args.filename}", values)


