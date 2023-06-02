import numpy as np
import os
import argparse

#  Beuatifully hardcoded script that scans the old Anderson data, where the gap is saved as the natural log
#  Saves an array with the 10-logarithm of the gap for each realization

parser = argparse.ArgumentParser()
parser.add_argument('--file_name', type=str)  # Name of file to save data to
args = parser.parse_args()

L_list = np.array([12, 24, 36, 48, 60, 70, 80])
W_list = np.array([0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
reals = 2000

values = np.zeros((len(L_list), len(W_list), reals))

for L_ind, L in np.ndenumerate(L_list):
    for W_ind, W in np.ndenumerate(W_list):
        print(f'L{L}, W{W}')
        path = f'data_anderson_old_only_spec/L{L}/W{W}'
        data_count = 0
        for entry in os.scandir(path):
            print(data_count)
            data = np.load(f'{path}/{entry.name}', allow_pickle=True)
            ln_gap = data[0]
            log10_gap = ln_gap / np.log(10)
            values[L_ind, W_ind, data_count] = log10_gap
            data_count += 1

np.save(f'npy_data/{args.file_name}', values)
