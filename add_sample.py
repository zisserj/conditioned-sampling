import argparse
import dd.cudd_add as _agd # type: ignore
import math
from time import process_time_ns
import numpy as np

import memray

np.set_printoptions(precision=3, suppress=True)
rng = np.random.default_rng()

from drdd_to_add import load_adds_from_drdd

_ms_str_from = lambda start_ns: f'{(process_time_ns()-start_ns)*1e-6:05.6f}ms'
_ms_str_any = lambda ns: f'{ns*1e-6:.6f}ms'


def make_sample_add(manager):
    two = manager.constant(2)
    ver_num = 4
    req_vars = math.ceil(math.log2(ver_num))
    x_var_names = [f'x{i}' for i in range(req_vars)]
    manager.declare(*x_var_names)
    y_var_names = [f'y{i}' for i in range(req_vars)]
    manager.declare(*y_var_names)
    z_var_names = [f'z{i}' for i in range(req_vars)]
    manager.declare(*z_var_names)
    x_vars = [manager.var(xi) for xi in x_var_names]
    y_vars = [manager.var(yi) for yi in y_var_names]
    opr_x = manager.apply('+', *x_vars) # doesnt scale for more vars but later problem
    opr_y = manager.apply('*', *y_vars)
    #opr_y = manager.apply('+', opr_y, two)
    opr_xy = manager.apply('\\/', opr_x, opr_y)
    return opr_xy

def _define_var_maps(ctx, trans):
    # find/define variables req
    x_var_names = []
    y_var_names = []
    for v in trans.support:
        if v.startswith('x'):
            x_var_names.append(v)
        elif v.startswith('y'):
            y_var_names.append(v)
        else:
            raise Exception("Unrecognised variable in g0")
    assert len(x_var_names) == len(y_var_names)
    ctx.var_length = len(x_var_names)
    x_var_names = sorted(x_var_names)
    y_var_names = sorted(y_var_names)
    z_var_names = [v.replace('x', 'z') for v in x_var_names]
    ctx.manager.declare(*z_var_names)

    map_mul = dict(zip(y_var_names, z_var_names)) | dict(zip(x_var_names, y_var_names))
    map_next_iter = dict(zip(z_var_names, y_var_names))
    return (map_mul, map_next_iter)
    
# trans must have all x1-xn, y1-yn vars
def compute_power_graphs(ctx, trans, n):
    map_mul, map_next_iter = _define_var_maps(ctx, trans)
    
    is_power_of_2 = (n & (n-1) == 0) and n != 0
    iters = int(np.log2(n))
    
    manager = ctx.manager
    gs = [trans]
    ts = []
    g_k = trans
    last_t = process_time_ns()
    for i in range(0, int(np.log2(n))):
        # t = g x g
        g_k_ = manager.let(map_mul, g_k)
        t_k = manager.apply('*', g_k, g_k_) # cuddGarbageCollect?
        ts.append(t_k)
        
        # if i <= iters or not is_power_of_2:
        # g = Ey in t
        g_k_pre = manager.exist(map_next_iter.values(), t_k)
        g_k = manager.let(map_next_iter, g_k_pre)
        gs.append(g_k)
        last_t = process_time_ns()
    return gs, ts

# bdd assignment as int repr
def _asgn_to_state(asgn, num_bits, vars=['x']):
    res = []
    for var in vars:
        x = ['1' if asgn.get(f'{var}{i}') else '0' for i in reversed(range(num_bits))]
        res.append(int(''.join(x), base=2))
    return res

# int repr as bdd assignment
def _state_to_asgn(state_ints, num_bits, vars):
    res = []
    for var, num in zip(vars, state_ints):
        bits = f'{num:0{num_bits}b}'[::-1]
        res.append({f'{var}{i}':True if b == '1' else False
                    for i, b in enumerate(bits)})
    return res

def _weighted_sample(opts_iter):
    coords, weights = zip(*opts_iter)
    weights = np.array(weights,dtype=float)
    weights /= weights.sum()
    return rng.choice(coords, axis=0, p=weights)

def _sample_add_conditioned(ctx, t, start, target, w):
    rename_map = {f'x{i}': f'z{i}' for i in range(ctx.var_length)}
    manager = ctx.manager
    target_rename = manager.let(rename_map, target)
    relevant_states = t & start & target_rename
    opts = list(manager.pick_iter(relevant_states, with_values=True))
    if len(opts) == 0:
        return "No matching traces"
    res = _weighted_sample(opts)
    res_ints = _asgn_to_state(res, ctx.var_length, 'xyz')
    w[0] = res_ints[0]
    w[len(w)//2] = res_ints[1]
    w[-1] = res_ints[2]

def _sample_seq_step(ctx, ti, lo, hi, w):
    manager = ctx.manager
    x_asgn, z_asgn = _state_to_asgn([w[lo],w[hi]], ctx.var_length, 'xz')
    start_bdd = manager.cube(x_asgn)
    target_bdd = manager.cube(z_asgn)
    relevant_states = ti & start_bdd & target_bdd
    opts = list(manager.pick_iter(relevant_states, with_values=True))
    res = _weighted_sample(opts)
    res_ints = _asgn_to_state(res, ctx.var_length, 'xyz')
    w[(lo+hi)//2] = res_ints[1]

def _draw_sample_fill(ctx, ts, t_idx, w, start, end):
    for i in range(t_idx, 0, -1):
        inc = np.power(2, i)
        for j in range(start, end, inc):
            _sample_seq_step(ctx, ts[i-1], j, j+inc, w)
    
def draw_sample_power(ctx, ts, length, init, target):
    w = [None]*(length+1)
    no_states = _sample_add_conditioned(ctx, ts[-1], init, target, w)
    if no_states:
        return no_states
    _draw_sample_fill(ctx, ts, int(np.log2(length))-1, w, 0, length)
    return w


def compute_forward_probs(ctx, gs, length, init):
    rename_map = {f'y{i}': f'x{i}' for i in range(ctx.var_length)}
    bin_rep = f'{length:b}'
    assert len(gs) >= len(bin_rep), f"Gs are missing for length {length}"
    
    prior_prob = init # would ideally divide by number of possible initial states but doesnt really matter
    steps_indices = []
    # forward compute
    for i, b in enumerate(reversed(bin_rep)): # lsb first
        if b == '1':
            # for all x in set prior_prop:={possible states to be at after prev gi steps starting at init} 
            # what is the probability of (having gotten to x) /\ (get to all y from x)
            yi_from_x = ctx.manager.apply('*', prior_prob, gs[i])
            
            # sum over x to get specific probability to be at state y after prev + cur gi
            prior_prob_y = ctx.manager.exist(rename_map.values(), yi_from_x)
            prior_prob = ctx.manager.let(rename_map, prior_prob_y)
            
            steps_indices.append((i, yi_from_x))
    return steps_indices

def draw_sample_generic(ctx, ts, length, target, mid_tranitions):
    rename_map = {f'x{i}': f'y{i}' for i in range(ctx.var_length)}
    
    w = [None]*(length+1)
    # backwards compute - given init and target, select middle nodes
    goal_idx = length
    
    # first iteration need to fill w
    steps_iter = reversed(mid_tranitions)
    i, trans_f = next(steps_iter)
    step_idx = goal_idx - (2**i)
    target_rename = ctx.manager.let(rename_map, target)
    
    opts = list(manager.pick_iter(trans_f & target_rename, with_values=True))
    if len(opts) == 0:
        return "No matching traces"
    res_pair =  _weighted_sample(opts)
    res_ints = _asgn_to_state(res_pair, ctx.var_length, 'xy')
    
    
    w[step_idx] = res_ints[0]
    w[goal_idx] = res_ints[1]
    _draw_sample_fill(ctx, ts, i, w, step_idx, goal_idx)
    goal_idx = step_idx
    
    for i, trans_f in steps_iter:
        # do endpoints sampling
        step_idx = goal_idx - (2**i)
        # goal state as assignment
        goal = {f'y{var[1:]}':val for var, val in res_pair.items() if var.startswith('x')}
        goal_bdd = ctx.manager.cube(goal)
        
        opts = list(manager.pick_iter(trans_f & goal_bdd, with_values=True))
        res_pair =  _weighted_sample(opts)
        res_ints = _asgn_to_state(res_pair, ctx.var_length, 'xy')
    
        w[step_idx] = res_ints[0]
        
        # "recursive" fill
        _draw_sample_fill(ctx, ts, i, w, step_idx, goal_idx)
        goal_idx = step_idx
    return w
 
def _state_to_og_vars(vars, w, intval):
    bits = f'{intval:0{w}b}'[::-1]
    idx= 0
    res = {}
    for name, num_bits in vars:
        res[name] = int(bits[idx:idx+num_bits], base=2)
        idx += num_bits
    return res

def generate_many_traces(ctx, gs, ts, length, init, target, save_traces=False, repeats=500, bypass=True):
    if not bypass and (length & (length-1) == 0) and length != 0: # https://stackoverflow.com/a/57025941
        draw = lambda: draw_sample_power(ctx, ts, length, init, target)
        #print(f"Property probability is {rel_mat.sum()/len(init)}")
    else:
        mid_transitions = compute_forward_probs(ctx, gs, length, init)
        draw = lambda: draw_sample_generic(ctx, ts, length, target, mid_transitions)
        # print(f"Property probability is {rel_mat.sum()/len(init)}")

    generated = []
    time_total = 0
    for _ in range(repeats):
        iter_start_time = process_time_ns()
        res = draw()
        time_total += process_time_ns() - iter_start_time
        if type(res) == str:
            print(res)
            return
        tr = tuple(res)
        if save_traces:
            generated.append(tr)
    ns_taken_avg = time_total / repeats
    print(f'Taken {_ms_str_any(ns_taken_avg)} per sample')
    if save_traces:
        return generated

def print_map():
    num_bits = len(transitions.support)//2
    g1_asgns = list(manager.pick_iter(gs[1], with_values=True))
    vars = [('s', 3), ('d', 3)]
    g1_dict = [_asgn_to_state(a, num_bits, "xy") for a, _ in g1_asgns]
    g1_v = [b for _, b in g1_asgns]
    g1_og = [{k: _state_to_og_vars(vars, 6, v) for k, v in entry.items()} for entry in g1_dict]
    print('\n'.join([f'{k}: {v}' for k, v in zip(g1_og, g1_v)]))

if __name__ == "__main__":
    
    parser = True
    # python add_sample.py dtmcs/dice/die.drdd 8 -repeats 10
    if parser:
        parser = argparse.ArgumentParser("Generates conditional samples of system via Algabraic Decision Diagrams.")
        parser.add_argument("fname", help="Model exported as drdd file by storm", type=str)
        parser.add_argument("length", help="Generated trace length (currently only supports powers of 2)", type=int)
        parser.add_argument("-repeats", help="Number of traces to generate", type=int, default=1000)
        parser.add_argument("-tlabel", help="Name of target label matching desired final states",
                            type=str, default='target')
        parser.add_argument('-output', help="File destination for generated traces", type=str, default='')
        parser.add_argument('--alg4', help="Use general alg for powers of two instead of baseline", action='store_true')
        parser.add_argument('--store', help="Store / try loading existing functions", action='store_true')
        args = parser.parse_args()
        filename = args.fname
        path_n = args.length
        repeats = args.repeats
        tlabel = args.tlabel
        store = args.store
        bypass = args.alg4
        output = args.output
    else:
        filename = "dtmcs/features/brp_N_32_MAX_3.drdd"
        path_n = 32
        repeats = 100
        tlabel = 'target'
        bypass = True
        store = False
        output = filename + '.out'
    print(f'Running parameters: fname={filename}, n={path_n},'+
          f' repeats={repeats}, label={tlabel}, store={store}, output={output if len(output) > 0 else False}')

    manager = _agd.ADD()
    manager.configure(max_growth=1.5)
    context = lambda: None # (required to assign attributes)
    context.manager = manager # type: ignore
    
    
    memlog = f'{filename}_{path_n}.bin'
    with memray.Tracker(destination=memray.FileDestination(memlog, overwrite=True), native_traces=True):
        parse_time = process_time_ns()
        model = load_adds_from_drdd(context.manager, filename)
        print(f'Finished parsing input: {_ms_str_from(parse_time)}.')
        init = model['initial']
        target = model[tlabel]
        assert len(target) > 0, "Target states missing"
        transitions = model['transitions']
        
        
        print(f"Number of variables per state: {len(transitions.support)//2}")
        print(f"Size of ADD: {transitions.dag_size} nodes")

        if store:
            raise NotImplementedError("ADD storage is not yet supported")
            # dirname = filename.replace('.drdd', '/')
            # gs, ts = load_and_store(dirname, transitions, path_n)
        else:
            precomp_time = process_time_ns()
            gs, ts = compute_power_graphs(context, transitions, path_n)
            print(f'Finished precomputing functions: {_ms_str_from(precomp_time)}.')
        
        save_traces = len(output) > 0
            
        res = generate_many_traces(context, gs, ts, path_n,
                    init, target, save_traces=save_traces,
                    repeats=repeats, bypass=bypass)
        
        if save_traces and res:
            with open(output, 'w+') as f:
                for label, states_bdd in model.items():
                    if label != 'transitions':
                        states_iter = context.manager.pick_iter(states_bdd)
                        states = [_asgn_to_state(s,context.var_length, 'x')[0]
                                    for s in states_iter]
                        #_asgn_to_state(res_pair, ctx.var_length, 'xy')
                        f.write(f'{label}: {str(states)}\n')
                f.write('---\n')
                f.write('\n'.join([str(r)[1:-1] for r in res]))
            print(f'{len(res)} traces written to {output}.')
        