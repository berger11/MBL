import numpy as np
import os
from statistics import mean

L_list = []
W_list = []
dir = "data"

for subdir, dirs, files in os.walk(dir):
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

for L_ind, L in enumerate(L_list):
    for W_ind, W in enumerate(W_list):

        path = f"./data/L{L}/W{W}"
        gap_list = []

        for entry in os.scandir(path):
            data = np.load(path+"/"+entry.name)
            gap_list.append(data[0])

        mean_gap = mean(gap_list)
        values[L_ind, W_ind] = mean_gap









#print(np.sort(L_list))
#print(np.sort(W_list))
