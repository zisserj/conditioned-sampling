import os
import time
import numpy as np
import scipy.sparse as sp
from drn_to_sparse import read_drn
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
# import matplot2tikz

np.set_printoptions(precision=5, suppress=True)
rng = np.random.default_rng()


_ms_str_from = lambda start_ns: f'{(time.perf_counter_ns()-start_ns)*1e-6:05.6f}ms'
_ms_str_any = lambda ns: f'{ns*1e-6:.6f}ms'

# make T[x,y,z] = G[x,y] * G[y,z]
def _compute_mid_step(g):
    wide = sp.block_diag(g, format='csr')
    mult = g @ wide
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
        ti = _compute_mid_step(gi)
        ts.append(ti)
        # print(ti.sum(axis=1))
    return gs, ts

def extend_power_mats(gs, ts, up_to):
    for i in range(len(gs)-1, up_to-1):
        gi = gs[i] @ gs[i]
        gs.append(gi)
    for i in range(len(ts), up_to):
        ti = _compute_mid_step(gs[i])
        ts.append(ti)

def plot_mats(dirname, gs, ts):
    fname = os.path.basename(dirname)
    n_func = len(gs)
    dim = gs[0].shape[0]
    msize = 500/(dim)
    cmap = mpl.colormaps['plasma']
    colors = cmap(np.linspace(0.8, 0, n_func))
    fig = plt.figure(figsize=(7,7))
    for i, g in enumerate(reversed(gs)):
        num_values = len(np.unique(g.data, sorted=False)) # type: ignore
        plt.spy(g, figure=fig, markersize=msize, marker='s',
                c=colors[i], label=f'$g_{n_func-i-1}$, {num_values} unique values')
    #plt.title(f'Matrix Rep of $g_i$ for "{fname}"')
    plt.xticks([0, dim-1])
    plt.yticks([0, dim-1])
    plt.legend(loc='lower left', markerscale=10/msize, reverse=True)
    #plt.tight_layout()
    plt.savefig(dirname+'_gs.png')
    # matplot2tikz.clean_figure()
    # matplot2tikz.save(f"{dirname}_gs.tex")
    print(f'exported {dirname}')
    plt.close()

def ts_sanity_test(ts, path_n, init, target):
        # P=? [F={path_n} "target"]
        actual_prob_test = _slice_csr_full(ts[int(np.log2(path_n))], init, target)
        print(actual_prob_test.sum())
        # should be the same as the storm property test
    
# M[i, j, k] = M'[i, (j*n)+k]
# arrays as input for inital/target states
# REALLY INEFFICIENT
def _slice_csr_full(mat, x, z):
    # every idx in z is a column in y
    d = mat.shape[0]
    if len(z) == 0:
        z = np.arange(d)
    cols = np.reshape(z, (-1, 1)) + np.arange(0, d**2, d)
    res = mat[x][:, cols.flatten()]
    newshape = res.reshape((len(x), -1, len(z)), order='F')
    # print(newshape.toarray())
    return newshape

# subsequent samples use simple indexing bc x, z arent arrays
def _slice_csr_col(mat, x, z):
    d = mat.shape[0]
    new_s = np.s_[x, z::d]
    return mat[new_s]

# pick a coordinate based on transition probability
def _weighted_idx_sample(mat):
    coords = np.array(mat.nonzero())
    weights = mat.data
    weights /= weights.sum() # normalize weights to 1
    res = rng.choice(coords, axis=1, p=weights)
    return res

# assumes initial states have the same probability of being chosen
def _sample_conditioned(ti, init, target, w, s=0, d=-1):
    d = len(w)+d if (d < 0) else d
    mid = (s+d)//2
    rel_mat = _slice_csr_full(ti, init, target)
    if rel_mat.sum() == 0:
        return "No matching traces"
    # per_init_idx = np.sum(rel_mat, axis=(0, 1))
    bounds_idx = _weighted_idx_sample(rel_mat)
    w[s] = init[bounds_idx[0]]
    w[mid] = bounds_idx[1]
    w[d] = target[bounds_idx[2]]


def _sample_seq_step(ti, lo, hi, w):
    mid = int(np.mean([lo,hi]))
    opts = _slice_csr_col(ti, w[lo], w[hi])
    asgn = _weighted_idx_sample(opts)
    w[mid] = asgn[0]
    # print(f'w[{mid}]={w[mid]}')

def _draw_sample_fill(ts, t_idx, w, start, end):
    for i in range(t_idx, 0, -1):
        inc = np.power(2, i)
        for j in range(start, end, inc):
            _sample_seq_step(ts[i-1], j, j+inc, w)

def draw_sample_simple(ts, length, init=[0], target=[]):
    w = np.full(length+1, -1, dtype=int)
    no_states = _sample_conditioned(ts[int(np.log2(length))-1], init, target, w)
    if no_states:
        return no_states
    _draw_sample_fill(ts, int(np.log2(length))-1, w, 0, length)
    return w

def compute_forward_probs(gs, length, init):
    n = gs[0].shape[0]
    bin_rep = f'{length:b}'
    assert len(gs) >= len(bin_rep), f"Gs are missing for length {length}"
    prior_prob = np.zeros(n)
    prior_prob[init] = 1/len(init) # uniform assumed for initial states
    # prior_prob = sp.csc_array(prior_prob)
    steps_indices = []
    # forward compute
    for i, b in enumerate(reversed(bin_rep)): # lsb first
        if b == '1':
            # for all x in set prior_prop:={possible states to be at after prev gi steps starting at init} 
            # what is the probability of (having gotten to x) /\ (get to all y from x)
            prior_diag = sp.diags_array(prior_prob, offsets=0)
            yi_from_x = prior_diag @ gs[i]
            
            # sum over x to get specific probability to be at state y after prev + cur gi
            prior_prob = yi_from_x.sum(axis=0)
            
            reachable = prior_prob.nonzero()
            steps_indices.append((i,
                                  sp.csr_array(yi_from_x)))
    return steps_indices

def draw_sample_generic(ts, length, target, mid_tranitions):
    w = np.full(length+1, -1, dtype=int)
    # backwards compute - given init and target, select middle nodes
    goal_idx = length
    
    # first iteration need to fill w
    steps_iter = reversed(mid_tranitions)
    i, trans_mat = next(steps_iter)
    step_idx = goal_idx - (2**i)
    end_pair =  _weighted_idx_sample(trans_mat[:,target])
    w[step_idx] = end_pair[0]
    w[goal_idx] = target[end_pair[1]]
    _draw_sample_fill(ts, i, w, step_idx, goal_idx)
    goal_idx = step_idx
    
    for i, trans_mat in steps_iter:
        # do endpoints sampling
        step_idx = goal_idx - (2**i)
        res_coord =  _weighted_idx_sample(trans_mat[:,w[goal_idx]])
        w[step_idx] = res_coord[0]
        # "recursive" fill
        _draw_sample_fill(ts, i, w, step_idx, goal_idx)
        #endpoint_sampled = [w[start_idx]]
        goal_idx = step_idx
    return w
    
def make_small_sample():
    dim = 4
    ts_data = np.array([[0,.3,0,0.7],[0,0.6,.4,0],[0,0,1,0],[0,1,0,0]])
    return sp.coo_matrix(ts_data, shape=(dim, dim), dtype=float).tocsr()


# alg works the same for counting number of traces, and is easier to see while testing
def make_small_sample_count():
    dim = 4
    ts_data = np.array([[0,1,0,1],[0,1,1,0],[0,0,1,0],[0,1,0,0]])
    return sp.coo_matrix(ts_data, shape=(dim, dim), dtype=float).tocsr()

def generate_many_traces(gs, ts, length, init, target, save_traces=False, repeats=500, bypass=True):
    init = np.array(init)
    target = np.array(target)
    if not bypass and (path_n & (path_n-1) == 0) and path_n != 0: # https://stackoverflow.com/a/57025941
        draw = lambda: draw_sample_simple(ts, length, init, target)
        rel_mat = _slice_csr_full(ts[-1], init, target)
        print(f"Property probability is {rel_mat.sum()/len(init)}")
    else:
        if len(gs) <= np.log2(path_n):
            extend_power_mats(gs, ts, len(gs)+1)
        endpoint_steps = compute_forward_probs(gs, length, init)
        draw = lambda: draw_sample_generic(ts, length, target, endpoint_steps)
        rel_mat = endpoint_steps[-1][1][:, target]
        print(f"Property probability is {rel_mat.sum()/len(init)}")
    generated = []
    time_total = 0
    
    for _ in range(repeats):
        iter_start_time = time.perf_counter_ns()
        res = draw()
        time_total += time.perf_counter_ns() - iter_start_time
        if type(res) == str:
            print(res)
            return
        tr = tuple(res.tolist())
        if save_traces:
            generated.append(tr)
    ns_taken_avg = time_total / repeats
    print(f'Taken {_ms_str_any(ns_taken_avg)} per sample')
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
        print(f'Finished precomputing remaining functions: {_ms_str_from(precomp_time)}.')
    else:
        precomp_time = time.perf_counter_ns()
        gs, ts = compute_power_mats(t0, length)
        print(f'Finished precomputing functions: {_ms_str_from(precomp_time)}.')
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
    # python sparse_mat_sample.py dtmcs/dice/die.drn 8 -repeats 10
    if parser:
        parser = argparse.ArgumentParser("Generates conditional samples of system via sparse matrices.")
        parser.add_argument("fname", help="Model exported as drn file by storm", type=str)
        parser.add_argument("length", help="Generated trace length", type=int)
        parser.add_argument("-repeats", help="Number of traces to generate", type=int, default=1000)
        parser.add_argument("-tlabel", help="Name of target label matching desired final states",
                            type=str, default='target')
        parser.add_argument('-output', help="File destination for generated traces", type=str, default='')
        parser.add_argument('--alg4', help="Use general alg for powers of two instead of baseline", action='store_true')
        parser.add_argument('--store', help="Store / try loading existing mats", action='store_true')
        args = parser.parse_args()
        filename = args.fname
        path_n = args.length
        repeats = args.repeats
        tlabel = args.tlabel
        bypass = args.alg4
        store = args.store
        output = args.output
    else:
        filename = "dtmcs/brp/brp_N_64_MAX_4.drn"
        path_n = 16
        repeats = 100
        tlabel = 'target'
        bypass = True
        store = False
        output = filename + '.out'
    print(f'Running parameters: fname={filename}, n={path_n}, repeats={repeats},'+
          f' label={tlabel}, store={store}, alg4={bypass}, output={output if len(output) > 0 else False}')
    parse_time = time.perf_counter_ns()
    model = read_drn(filename)
    print(f'Finished parsing input: {_ms_str_from(parse_time)}.')
    init = model['init']
    assert tlabel in model, f"Target label '{tlabel}' missing"
    target = model[tlabel]
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
        print(f'Finished precomputing functions: {_ms_str_from(precomp_time)}.')
    
    
    save_traces = len(output) > 0
    # plot_mats(filename.replace('.drn', ''), gs, ts)
    # quit(0)
    res = generate_many_traces(gs, ts, path_n, init,
                target, repeats=repeats, save_traces=save_traces, bypass=bypass)
    
    if save_traces and res:
        with open(output, 'w+') as f:
            for label, states in model.items():
                if label != 'trans':
                    f.write(f'{label}: {str(states)}\n')
            f.write('---\n')
            f.write('\n'.join([str(r)[1:-1] for r in res]))
        print(f'{len(res)} traces written to {output}.')

