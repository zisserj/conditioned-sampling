import numpy as np
from itertools import product, chain
import pickle
import argparse
import scipy.sparse as sp

parser = argparse.ArgumentParser("Process transitions BDD into adjacency matrix.")
# parser.add_argument("fname", help="Target transitions (BDD) file.", type=str)
# args = parser.parse_args()
# fname = args.fname
fname = "smv_examples/semaphore_trans.bdd"

with open(fname) as f:
    content = [s.strip() for s in f.readlines()]

sep = content.index('')
vars = np.array(content[:sep], dtype=np.str_)
table = content[sep+1:]

vars_next = np.strings.startswith(vars, 'next(')



table[0] = 'TRUE'
table[1] = 'FALSE'
for i in range(2, len(table)):
    table[i] = table[i].strip().split()
    table[i] = [int(j) for j in table[i]]

T_asgn = [None for _ in table]
T_asgn[0] = []
F_asgn = [None for _ in table]
F_asgn[1] = []

def add_to_assignment(var, asgn, dicts):
    if dicts is None:
        return None
    return ((var, asgn), dicts)

def extend_if_exists(cur, new_list):
    if cur and new_list:
        cur.extend(new_list)
        return cur
    elif new_list:
        return new_list
    return cur

for i, row in enumerate(table[2:], start=2):
    var = row[0]
    new_true = []
    new_false = []
    new_true.append(add_to_assignment(var, True, T_asgn[row[2]]))
    new_false.append(add_to_assignment(var, True, F_asgn[row[2]]))
    # low edge complemented
    if row[1] == 1:
        new_true.append(add_to_assignment(var, False, F_asgn[row[3]]))
        new_false.append(add_to_assignment(var, False, T_asgn[row[3]]))
    else:
        new_true.append(add_to_assignment(var, False, T_asgn[row[3]]))
        new_false.append(add_to_assignment(var, False, F_asgn[row[3]]))
    T_asgn[i] = [a for a in new_true if a != None]
    F_asgn[i] = [a for a in new_false if a != None]


def as_string(assignment):
    s = ''
    for i in range(len(vars)):
        if i in assignment:
            c = '■' if assignment[i] else '□'
        else:
            c = '▣'
        s += c
    return s

true_as = [as_string(a) for a in T_asgn[-1]]
false_as = [as_string(a) for a in F_asgn[-1]]
true_as.sort()
false_as.sort()

def flatten_assignments(t):
    if len(t[1]) == 0:
        asgn = [None]*len(vars)
        asgn[t[0][0]] = t[0][1]
        return [asgn]
    asgns = chain.from_iterable([flatten_assignments(r) for r in t[1]])
    res = []
    for a in asgns:
        a[t[0][0]] = t[0][1]
        res.append(a)
    return res

# https://stackoverflow.com/a/64391584
def binarr_to_deci(binary):
    return sum(val*(2**idx) for idx, val in enumerate(reversed(binary)))

def bin_to_tuple(a):
    return binarr_to_deci(a[::2]), binarr_to_deci(a[1::2])

def rec_asgn_as_arr(t):
    base_arr = np.array(t)
    base_arr = np.array(base_arr)
    opts = base_arr == None
    if sum(opts) == 0:
        return [bin_to_tuple(base_arr)]
    res = []
    for c in product([False,True], repeat=sum(opts)):
        base_arr[opts] = c
        res.append(bin_to_tuple(base_arr))
    return res

def asgn_as_arr(dict):
    base_arr = [None]*(len(vars))
    for k, v in dict.items():
        base_arr[k] = '1' if v else '0'
    to_tuple = lambda a: (int(''.join(a[::2]), base=2), int(''.join(a[1::2]), base=2))
    base_arr = np.array(base_arr)
    opts = base_arr == None
    if sum(opts) == 0:
        return [to_tuple(base_arr)]
    res = []
    for c in product('01', repeat=sum(opts)):
        base_arr[opts] = c
        res.append(to_tuple(base_arr))
    return res


def dict_to_mat(true_asgns):
    dim = 2**(len(vars)//2)
    mat = sp.dok_array((dim, dim),dtype=int)
    flatened = list(chain.from_iterable([flatten_assignments(t) for t in true_asgns]))
    for sat_asgn in flatened:
        for (i, j) in rec_asgn_as_arr(sat_asgn):
            mat[i,j] = 1
    return mat

# print([[0]*4 for _ in range(len(vars)//2)])
res = dict_to_mat(T_asgn[-1])
print(res.toarray())

with open(fname.replace('.bdd', '.npz'), 'wb+') as f:
    pickle.dump(res, f)