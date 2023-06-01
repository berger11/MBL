import numpy as np
import os
from statistics import mean
import argparse
from functions.flint_help_functions import string_to_acb
from flint import arb, ctx
# Scans data from a given directory and computes the 10-log of the gap, and saves each realization

ctx.dps = 300

parser = argparse.ArgumentParser()
parser.add_argument('--sim_name', type=str)  # Name of simulation to scan data from
parser.add_argument('--file_name', type=str)  # Name of file to save data to
args = parser.parse_args()
L_list = []
W_list = []

for subdir, dirs, files in os.walk(f'data/{args.sim_name}_data'):
    print(subdir)
    subdir_split = subdir.split('/')
    last = subdir_split[-1]
    if last[0] == 'L':
        number = last[1:]
        if int(number) not in L_list:
            L_list.append(int(number))
    if last[0] == 'W':
        number = last[1:]
        if float(number) not in W_list:
            W_list.append(float(number))

values = np.zeros((len(L_list), len(W_list)), dtype=object)
L_list = np.sort(L_list)
W_list = np.sort(W_list)
print(L_list)
print(W_list)
for L_ind, L in enumerate(L_list):
    for W_ind, W in enumerate(W_list):
        print(f'L{L} W{W}')
        path = f'data/{args.sim_name}_data/L{L}/W{W}'
        data_count = 0
        real_list = []
        for entry in os.scandir(path):
            data = np.load(f'{path}/{entry.name}', allow_pickle=True)
            config = data[-1]
            vectors = config['vectors']
            eigen = data[0]
            if vectors:
                spectrum = eigen[0]
            else:
                spectrum = eigen
            try:
                assert L == config['length']
                assert W == config['disorder']
            except AssertionError:
                print('WARNING: DIRECTORY STRUCTURE AND CONFIG FILES NOT IN AGREEMENT')
            for i in range(4*L):
                spectrum[i] = string_to_acb(spectrum[i])  # Converts the eigenvalues from strings back to acb

            gap = np.abs(np.real(spectrum[0]))
            for j in range(len(spectrum)):
                if 0 < np.real(spectrum[j]) < gap:
                    gap = np.real(spectrum[j])
            gap = 2 * gap
            log_gap = float(gap.log() / arb(10).log())

            print('log gap', log_gap)
            real_list.append(log_gap)
            data_count += 1
        try:
            assert data_count == config['reals']
        except AssertionError:
            print('WARNING: NUMBER OF FILES NOT EQUAL TO NUMBER OF REALIZATIONS IN CONFIG FILE')

        assert values[L_ind, W_ind] == 0
        real_list = np.array(real_list)
        values[L_ind, W_ind] = real_list

np.save(f'npy_data/{args.file_name}', values)


