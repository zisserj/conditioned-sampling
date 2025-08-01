import re
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

pat = r"--- \w+\/([a-zA-Z]+)_([\S_]+).drn - (\d+) ---\
[\w ]+input: ([\d.]+)ms.\
[\w ]+states: ([\d.]+)\
[\w ]+transitions: ([\d.]+)\
[\w ]+functions: ([\d.]+)ms.\
[#\d\-\s]+Taken ([\d.]+)ms per sample"

# model, params, length, parse, #states, #transitions, precompute, avg/sample

fname = "/home/jules/conditioned_sampling/dtmcs/logs/mat_sampling-all-5863481.out"
with open(fname) as f:
    content = f.read()

res = [] # length, states, transitinons, precomp, sample
for match in re.finditer(pat, content):
    print(match.groups())
    #res.append([match.group(i) for i in [1, 2, 3, 4, 5, 6, 7, 8]])

# with open('dtmcs/brp_temp.csv', 'w') as f:
#     f.write('length, states, trans, precomp, sample\n')
#     f.writelines([(', '.join([p for p in line]) + '\n') for line in res])
