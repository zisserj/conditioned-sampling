import re
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# pat = r"--- \w+\/([a-zA-Z]+)_([\S_]+).drn - (\d+) ---\
# [\w ]+input: ([\d.]+)ms.\
# [\w ]+states: ([\d.]+)\
# [\w ]+transitions: ([\d.]+)\
# [\w ]+functions: ([\d.]+)ms.\
# [#\d\-\s]*Taken ([\d.]+)ms per sample"
pat = r"--- \w+\/([a-zA-Z]+)_?([\S_]+).drn - (\d+) ---\
[\s\S]*?\
[\w ]+input: ([\d.]+)ms.\
[\w ]+states: ([\d.]+)\
[\w ]+transitions: ([\d.]+)\
[\s\S]*?\
?[\w ]+functions: ([\d.]+)ms.\
[\w ]+ mats: G\(\d+\), T\((\d+)\)\
[\w ]+ probability is ([\d.+-e]+)\
[#\d\-\s]*Taken ([\d.]+)ms per sample"



content = ""
head = "dtmcs/logs/"
fnames = ["mat_sampling-brp-5994883.out",
          "mat_sampling-herman9-5994914.out",
          "mat_sampling-crowds-5994925.out",
          "mat_sampling-leader-5994880.out",
          "mat_sampling-herman-5994866.out"]
for fname in fnames:
    with open(head+fname) as f:
        content += f.read()

res = [] # length, states, transitinons, precomp, sample
for match in re.finditer(pat, content):
    print(match.groups())
    res.append(match.groups())
    #res.append([match.group(i) for i in [2, 3, 4, 5, 6, 7, 8]])

# model, params, length, parse, #states, #transitions,
# precompute, written mats, prob, avg/sample

with open('dtmcs/new_mat_timing.csv', 'w') as f:
    f.write('model, params, length, parsetime, states, trans, precomptime, newmat, prob, sampletime\n')
    f.writelines([(', '.join([p for p in line]) + '\n') for line in res])
