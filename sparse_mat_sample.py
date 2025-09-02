import os
import time
import numpy as np
import scipy.sparse as sp
import itertools
from add_sample import draw_sample
from drn_to_sparse import read_drn
import argparse

np.set_printoptions(precision=2, suppress=True)
rng = np.random.default_rng()


ms_str_from = lambda start_ns: f'{(time.perf_counter_ns()-start_ns)*1e-6:05.6f}ms'
ms_str_any = lambda ns: f'{ns*1e-6:.6f}ms'


# make T[x,y,z] = G[x,y] * G[y,z]
def compute_mid_step(g):
    wide = sp.block_diag(g)
    mult = g.tocoo() @ wide.tocoo()
    # print(mult.data.nbytes + mult.indptr.nbytes + mult.indices.nbytes)
    return mult

def compute_power_mats(trans, length):
    gs = [trans]
    gi = trans
    for i in range(1, int(np.log2(length))):
        gi = gi @ gi
        gs.append(gi)
    ts = []
    for gi in gs:
        ti = compute_mid_step(gi)
        ts.append(ti)
        # print(ti.sum(axis=1))
    return gs, ts

def extend_power_mats(gs, ts, up_to):
    for i in range(len(gs)-1, up_to-1):
        gi = gs[i] @ gs[i]
        gs.append(gi)
    for i in range(len(ts), up_to):
        ti = compute_mid_step(gs[i])
        ts.append(ti)

def ts_sanity_test(ts, path_n, init, target):
        # P=? [F={path_n} "target"]
        actual_prob_test = slice_csr_full(ts[int(np.log2(path_n))], init, target)
        print(actual_prob_test.sum())
        # should be the same as the storm property test
    
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
def sample_conditioned(ti, init, target, w, s=0, d=-1):
    d = len(w)+d if (d < 0) else d
    mid = (s+d)//2
    rel_mat = slice_csr_full(ti, init, target)
    if np.isclose(rel_mat.max(), 0):
        return "No matching traces"
    # per_init_idx = np.sum(rel_mat, axis=(0, 1))
    bounds_idx = weighted_idx_sample(rel_mat)
    w[s] = init[bounds_idx[0]]
    w[mid] = bounds_idx[1]
    w[d] = target[bounds_idx[2]]


def sample_seq_step(ti, lo, hi, w):
    mid = int(np.mean([lo,hi]))
    opts = slice_csr_col(ti, w[lo], w[hi])
    asgn = weighted_idx_sample(opts)
    w[mid] = asgn[0]
    # print(f'w[{mid}]={w[mid]}')

def draw_sample_fill(ts, t_idx, w, start, end):
    for i in range(t_idx, 0, -1):
        inc = np.power(2, i)
        for j in range(start, end, inc):
            sample_seq_step(ts[i-1], j, j + inc, w)

def draw_sample_simple(ts, length, init=[0], target=[]):
    w = np.full(length+1, -1, dtype=int)
    no_states = sample_conditioned(ts[-1], init, target, w)
    if no_states:
        return no_states
    draw_sample_fill(ts, int(np.log2(length))-1, w, 0, length)
    return w

def compute_nonpower_indices(gs, length, init):
    bin_rep = f'{length:b}'
    assert len(gs) >= len(bin_rep), f"Gs are missing for length {length}"
    reachable = init
    prior_prob = np.ones(gs[0].shape[0])[init]/len(init) # uniform assumed for initial states
    steps_indices = []
    # forward compute
    for i, b in enumerate(reversed(bin_rep)): # lsb first
        if b == '1':
            # for all x in set prior_init:={possible states to be at after prev gi steps starting at init} 
            # what is the probability of (having gotten to x) /\ (get to all y from x)
            yi_from_x = gs[i][reachable].multiply(prior_prob.reshape(-1,1))
            # sum over x to get specific probability to be at state y after prev + cur gi
            marginal_yi = yi_from_x.sum(axis=0)
            marginal_yi = marginal_yi/marginal_yi.sum() # normalize to get probabilities (its all relative anyway)
            
            reachable = marginal_yi.nonzero()
            prior_prob = marginal_yi[reachable]
            steps_indices.append((i,
                                  sp.csr_array(marginal_yi)))
    return steps_indices

def draw_sample_nonpower(gs, ts, length, init, target, indices):
    w = np.full(length+1, -1, dtype=int)
    # backwards compute - given init and target, select middle nodes
    steps_iter = reversed(indices)
    g_idx, rel_states = next(steps_iter)
    try:
        endpoint_sampled = target[weighted_idx_sample(rel_states[target])]
    except:
        return "No matching traces"
    target_idx = length
    start_idx = target_idx - (2**(g_idx))
    for prev_g_idx, prev_state in steps_iter:
        # do highest order sampling
        sample_conditioned(ts[g_idx-1], prev_state.nonzero()[0],
                        endpoint_sampled, w, s=start_idx, d=target_idx)
        # "recursive" fill
        draw_sample_fill(ts, g_idx-1, w, start_idx, target_idx)
        g_idx = prev_g_idx
        endpoint_sampled = [w[start_idx]]
        target_idx = start_idx
        start_idx -= 2**g_idx
    if g_idx > 0: # last step is at least 2
        sample_conditioned(ts[g_idx], init,
                        endpoint_sampled, w, s=start_idx, d=target_idx)
        # "recursive" fill
        draw_sample_fill(ts, g_idx, w, start_idx, target_idx)
    else: # last step is exactly 1
        opts = sp.coo_array(gs[0][init, endpoint_sampled])
        init_idx = weighted_idx_sample(opts)
        w[0] = init[init_idx][0] # need to include g0 in args?
    return w


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

def generate_many_traces(gs, ts, length, init, target, save_traces=False, repeats=500):
    if np.log2(path_n) == np.floor(np.log2(path_n)):
        draw = lambda: draw_sample_simple(ts, length, init, target)
    else:
        if len(gs) < np.log2(path_n):
            extend_power_mats(gs, ts, len(gs)+1)
        g_steps = compute_nonpower_indices(gs, length, init)
        draw = lambda: draw_sample_nonpower(gs, ts, length, init, target, g_steps)
    generated = []
    time_total = 0
    rel_mat = slice_csr_full(ts[-1], init, target)
    print(f"Property probability is {rel_mat.sum()/len(init)}")
    for _ in range(repeats):
        iter_start_time = time.perf_counter_ns()
        res = draw()
        if type(res) == str:
            print(res)
            return
        tr = tuple(res.tolist())
        time_total += time.perf_counter_ns() - iter_start_time
        if save_traces:
            generated.append(tr)
    ns_taken_avg = time_total / repeats
    print(f'Taken {ms_str_any(ns_taken_avg)} per sample')
    if save_traces:
        return generated
    
def load_and_store(dirname, t0, length):
    os.makedirs(dirname, exist_ok=True)
    num_mats = int(np.log2(length))
    gs, ts = [], []
    mat_fname = dirname + '{}{}.npz'
    for i in range(num_mats):
        if os.path.exists(mat_fname.format('G', i)):
            gs.append(sp.load_npz(mat_fname.format('G', i)))
        else:
            break
        if os.path.exists(mat_fname.format('T', i)):
            ts.append(sp.load_npz(mat_fname.format('T', i)))
        else:
            break
    exist_gs = len(gs)
    exist_ts = len(ts)
    if exist_ts == num_mats:
        print(f'Found all required mats.')
        return gs, ts
    elif exist_gs + exist_ts > 0:
        print(f'Found prior mats: G({exist_gs-1}), T({exist_ts-1})')
        precomp_time = time.perf_counter_ns()
        extend_power_mats(gs, ts, num_mats)
        print(f'Finished precomputing remaining functions: {ms_str_from(precomp_time)}.')
    else:
        precomp_time = time.perf_counter_ns()
        gs, ts = compute_power_mats(t0, length)
        print(f'Finished precomputing functions: {ms_str_from(precomp_time)}.')
    for i in range(exist_gs, len(gs)):
        mat_G = dirname + 'G{}.npz'
        sp.save_npz(mat_G.format(i), gs[i])
    for i in range(exist_ts, len(ts)):
        mat_T = dirname + 'T{}.npz'
        sp.save_npz(mat_T.format(i), ts[i])
    print(f'Stored generated mats: G({len(gs)-1}), T({len(ts)-1})')
    return gs, ts


if __name__ == "__main__":
    parser = True
    # python sparse_mat_sample.py dtmcs/die.drn 8 -repeats 10
    if parser:
        parser = argparse.ArgumentParser("Generates conditional samples of system via sparse matrices.")
        parser.add_argument("fname", help="Model exported as drn file by storm", type=str)
        parser.add_argument("length", help="Generated trace length (currently only supports powers of 2)", type=int)
        parser.add_argument("-repeats", help="Number of traces to generate", type=int, default=1000)
        parser.add_argument("-tlabel", help="Name of target label matching desired final states",
                            type=str, default='target')
        parser.add_argument('--store', help="Store / try loading existing mats", action='store_true')
        args = parser.parse_args()
        filename = args.fname
        path_n = args.length
        repeats = args.repeats
        tlabel = args.tlabel
        store = args.store
    else:
        filename = "dtmcs/die.drn"
        path_n = 8
        repeats = 100
        tlabel = 'target'
        store = False
    print(f'Running parameters: fname={filename}, n={path_n}, repeats={repeats}, label={tlabel}, store={store}')
    parse_time = time.perf_counter_ns()
    model = read_drn(filename, target_label=tlabel)
    print(f'Finished parsing input: {ms_str_from(parse_time)}.')
    init = model['init']
    target = model['target']
    assert len(target) > 0, "Target states missing"
    transitions = model['trans'].tocsr()
    
    print(f"Number of states: {transitions.shape[0]}")
    print(f"Number of transitions: {transitions.nnz}")
    
    if store:
        dirname = filename.replace('.drn', '/')
        gs, ts = load_and_store(dirname, transitions, path_n)
    else:
        precomp_time = time.perf_counter_ns()
        gs, ts = compute_power_mats(transitions, path_n)
        print(f'Finished precomputing functions: {ms_str_from(precomp_time)}.')
    
    
    # print(f'Finished drawing 1 sample: {ms_from(precomp_time)} from parse')
    res = generate_many_traces(gs, ts, path_n, init, target, repeats=repeats)
     
    
