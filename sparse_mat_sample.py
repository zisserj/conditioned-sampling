# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse as sp
import itertools
np.set_printoptions(precision=2, suppress=True)
rng = np.random.default_rng()
from drn_to_sparse import read_drn
# import argparse

# parser = argparse.ArgumentParser("Precomutes and uniformly samples TS via matrix representaation.")
# parser.add_argument("fname", help="Target mat pickle.", type=str)
# args = parser.parse_args()
# fname = args.fname


# make T[x,y,z] = G[x,y] * G[y,z]
def compute_mid_step(g):
    wide = sp.block_diag(g)
    mult = g @ wide
    return mult

def generate_power_mats(transition, length):
    gs = [transition]
    gi = transition
    for i in range(1, int(np.log2(length))):
        gi = gi @ gi
        gs.append(gi)
    #     print(f"i = {i}")

    ts = []
    for gi in gs:
        ti = compute_mid_step(gi)
        ts.append(ti)
        # print(ti.sum(axis=1))
    return gs, ts

# M[i, j, k] = M'[i, (j*n)+k]
# arrays as input for inital/target states
def slice_csr_full(mat, x, z):
    # every idx in z is a column in y
    d = mat.shape[0]
    if len(z) == 0:
        z = np.arange(d)
    cols = np.reshape(z, (-1, 1)) + np.arange(0, d**2, d)
    res = mat[x][:, cols.flatten()]
    newshape = res.reshape((len(x), -1, len(z)), order='F')
    # print(newshape.toarray())
    return newshape

# it gave me so much grief...
def slice_csr_bad(mat, x, z):
    # every idx in z is a column in y
    d = mat.shape[0]
    if len(z) == 0:
        z = np.arange(d)
    cols = (np.reshape(z, (-1, 1)) * d) + np.arange(d)
    res = mat[x][:, cols.flatten()]
    newshape = res.reshape((len(x), len(z), -1))
    print(newshape.toarray())
    return newshape

# subsequent samples use simple indexing bc x, z arent arrays
def slice_csr_col(mat, x, z):
    d = mat.shape[0]
    new_s = np.s_[x, z::d]
    return mat[new_s]

def test_slice_csr(actual):
    goal = np.array([[ 3472,  3808,  4144],
                    [10788, 12152, 13516],
                    [19448, 22032, 24616],
                    [29452, 33448, 37444]])
    for i in range(goal.shape[0]):
        for j in range(goal.shape[1]):
            for k in range(goal.shape[2]):
                is_eq = goal[i,j,k] == actual[i,j,k]
                if not is_eq:
                    print(i,j,j)


# pick a coordinate based on transition probability
def weighted_idx_sample(mat):
    coords = np.array(mat.nonzero())
    weights = mat.data
    weights /= weights.sum() # normalize weights to 1
    res = rng.choice(coords, axis=1, p=weights)
    return res


# assumes initial states have the same probability of being chosen
def sample_conditioned(ti, init, target, w):
    mid = (len(w))//2
    rel_mat = slice_csr_full(ti, init, target)
    # per_init_idx = np.sum(rel_mat, axis=(0, 1))
    bounds_idx = weighted_idx_sample(rel_mat)
    w[0] = init[bounds_idx[0]]
    w[mid] = bounds_idx[1]
    w[-1] = target[bounds_idx[2]]

def sample_seq_step(ti, lo, hi, w):
    mid = int(np.mean([lo,hi]))
    opts = slice_csr_col(ti, w[lo], w[hi])
    asgn = weighted_idx_sample(opts)
    w[mid] = asgn[0]
    # print(f'w[{mid}]={w[mid]}')

def draw_sample(ts, length, init=[0], target=[]):
    w = np.full(length+1, -1, dtype=int)
    sample_conditioned(ts[-1], init, target, w)
    for i in range(int(np.log2(path_n))-1, 0, -1):
        inc = np.power(2, i)
        for j in range(0, path_n, inc):
            sample_seq_step(ts[i-1], j, j + inc, w)
    return w

def draw_sample_init(gs, length, init, target):
    bin_rep = f'{length:b}'
    assert len(gs) >= len(bin_rep), f"Gs are missing for length {length}"
    reachable = [init]
    prior_prob = [np.ones(dim)[init]/len(init)]
    for i, b in enumerate(reversed(bin_rep)): # lsb first
        if b == '1':
            # for all x in set prior_init:={possible states to be at after prev gi steps starting at init} 
            # what is the probability of (having gotten to x) /\ (get to all y from x)
            yi_from_x = gs[i][reachable[-1]].multiply(prior_prob[-1].reshape(-1,1))
            # sum over x to get specific probability to be at state y after prev + cur gi
            marginal_yi = yi_from_x.sum(axis=0)
            #res = res/res.sum() # normalize to get probabilities (its all relative anyway)
            reachable.append(marginal_yi.nonzero())
            prior_prob.append(marginal_yi[reachable[-1]])
    print(marginal_yi[target]) # type: ignore

def make_small_sample():
    dim = 4
    ts_data = [0,.3,0,0.7,0,0.6,.4,0,0,0,1,0,0,1,0,0]
    ts_t0 = {k: v for k, v in zip(itertools.product(range(4), repeat=2), ts_data) if v != 0}
    row = [i for (i, j) in ts_t0.keys()]
    col = [j for (i, j) in ts_t0.keys()]
    vals = list(ts_t0.values())
    return sp.coo_array((vals, (row, col)), shape=(dim, dim), dtype=float).tocsr()


# alg works the same for counting number of traces, and is easier to see while testing
def make_small_sample_count():
    dim = 4
    ts_data = [0,1,0,1,0,1,1,0,0,0,1,0,0,2,0,0]
    ts_t0 = {k: v for k, v in zip(itertools.product(range(4), repeat=2), ts_data) if v != 0}

    row = [i for (i, j) in ts_t0.keys()]
    col = [j for (i, j) in ts_t0.keys()]
    vals = list(ts_t0.values())
    return sp.coo_array((vals, (row, col)), shape=(dim, dim), dtype=float).tocsr()

def generate_many_traces(ts, init, target, repeats=1000):
    results = {}
    for _ in range(repeats):
        tr = tuple(draw_sample(ts, path_n, init, target))
        if tr not in results:
            results[tr] = 1
        else:
            results[tr] += 1
    legible = '\n'.join([','.join([str(i) for i in k]) + f' - {v}' for k, v in results.items()])
    print(legible)


if __name__ == "__main__":
    filename = "/home/jules/storm_sampler/storm-project-starter-cpp/sparse_model.drn"

    model = read_drn(filename)
    init = model['init']
    goal = model['target']
    transitions = model['trans']
    # transitions = make_small_sample()
    
    # goal = [1,2]
    
    path_n = 8
    dim = transitions.shape[0] # type: ignore
    
    gs, ts = generate_power_mats(transitions, path_n)
    trace = draw_sample(ts, path_n, init, goal)
    # print(trace)
    generate_many_traces(ts, init, goal)
    
     
    #res = draw_sample_init(gs, 6, [0,1], [0,1,2])
    
    
    
