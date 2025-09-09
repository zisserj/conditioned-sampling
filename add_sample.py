import argparse
import dd.cudd_add as _agd # type: ignore
import math
from time import perf_counter_ns
import numpy as np

np.set_printoptions(precision=2, suppress=True)
rng = np.random.default_rng()

from drdd_to_add import load_adds_from_drdd
from bdd_prob_sample import state_to_og_vars # type: ignore

ms_str_from = lambda start_ns: f'{(perf_counter_ns()-start_ns)*1e-6:05.6f}ms'
ms_str_any = lambda ns: f'{ns*1e-6:.6f}ms'


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

def define_var_maps(ctx, trans):
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
def compute_power_graphs(ctx, trans, length):
    map_mul, map_next_iter = define_var_maps(ctx, trans)
    
    manager = ctx.manager
    gs = [trans]
    ts = []
    g_k = trans
    last_t = perf_counter_ns()
    for i in range(0, int(np.log2(length))):
        # t = g x g
        g_k_ = manager.let(map_mul, g_k)
        t_k = manager.apply('*', g_k, g_k_) # cuddGarbageCollect?
        ts.append(t_k)
        
        if i < int(np.log2(length))-1:
            # g = Ey in t
            g_k_pre = manager.exist(map_next_iter.values(), t_k)
            g_k = manager.let(map_next_iter, g_k_pre)
            gs.append(g_k)
        #print(manager.statistics())
        # print(f'Finished iteration {i}: {(perf_counter_ns()-last_t)*1e-9}')
        last_t = perf_counter_ns()
    return gs, ts

def asgn_to_state(asgn, num_bits, vars=['x']):
    res = []
    for var in vars:
        x = ['1' if asgn.get(f'{var}{i}') else '0' for i in reversed(range(num_bits))]
        res.append(int(''.join(x), base=2))
    return res

def state_to_asgn(state_ints, num_bits, vars):
    res = []
    for var, num in zip(vars, state_ints):
        bits = f'{num:0{num_bits}b}'[::-1]
        res.append({f'{var}{i}':True if b == '1' else False
                    for i, b in enumerate(bits)})
    return res

def weighted_sample(opts_iter):
    coords, weights = zip(*opts_iter)
    #coords = np.array(coords)
    weights = np.array(weights,dtype=float)
    weights /= weights.sum()
    return rng.choice(coords, axis=0, p=weights)

def sample_add_conditioned(ctx, t, start, target, w):
    rename_map = {f'x{i}': f'z{i}' for i in range(ctx.var_length)}
    manager = ctx.manager
    target_rename = manager.let(rename_map, target)
    rel_states = t & start & target_rename
    opts = list(manager.pick_iter(rel_states, with_values=True))
    if len(opts) == 0:
        return "No matching traces"
    res = weighted_sample(opts)
    res_ints = asgn_to_state(res, ctx.var_length, 'xyz')
    w[0] = res_ints[0]
    w[len(w)//2] = res_ints[1]
    w[-1] = res_ints[2]

def sample_bdd_seq(ctx, t, x_idx, z_idx, w):
    manager = ctx.manager
    x_asgn, z_asgn = state_to_asgn([w[x_idx],w[z_idx]], ctx.var_length, 'xz')
    start_bdd = manager.cube(x_asgn)
    target_bdd = manager.cube(z_asgn)
    rel_states = t & start_bdd & target_bdd
    opts = list(manager.pick_iter(rel_states, with_values=True))
    res = weighted_sample(opts)
    res_ints = asgn_to_state(res, ctx.var_length, 'xyz')
    w[(x_idx+z_idx)//2] = res_ints[1]

def draw_sample(ctx, ts, length, init, target):
    w = [None]*(length+1)
    no_states = sample_add_conditioned(ctx, ts[-1], init, target, w)
    if no_states:
        return no_states
    for i in range(int(np.log2(length))-1, 0, -1):
        inc = np.power(2, i)
        for j in range(0, length, inc):
            sample_bdd_seq(context, ts[i-1], j, j + inc, w)
    return w
    
def state_to_og_vars(vars, w, intval):
    bits = f'{intval:0{w}b}'[::-1]
    idx= 0
    res = {}
    for name, num_bits in vars:
        res[name] = int(bits[idx:idx+num_bits], base=2)
        idx += num_bits
    return res

def generate_many_traces(ctx, ts, length, init, target, save_traces=False, repeats=500):
    generated = []
    time_total = 0
    # todo: print probability of property
    #print(f"Property probability is {rel_mat.sum()/len(init)}")
    for _ in range(repeats):
        iter_start_time = perf_counter_ns()
        res = draw_sample(ctx, ts, length, init, target)
        if type(res) == str:
            print(res)
            return
        tr = tuple(res)
        time_total += perf_counter_ns() - iter_start_time
        if save_traces:
            generated.append(tr)
    ns_taken_avg = time_total / repeats
    print(f'Taken {ms_str_any(ns_taken_avg)} per sample')
    if save_traces:
        return generated

def print_graph():
    num_bits = len(transitions.support)//2
    g1_asgns = list(manager.pick_iter(gs[1], with_values=True))
    vars = [('s', 3), ('d', 3)]
    g1_dict = [asgn_to_state(a, num_bits, "xy") for a, _ in g1_asgns]
    g1_v = [b for _, b in g1_asgns]
    g1_og = [{k: state_to_og_vars(vars, 6, v) for k, v in entry.items()} for entry in g1_dict]
    print('\n'.join([f'{k}: {v}' for k, v in zip(g1_og, g1_v)]))

if __name__ == "__main__":
    
    parser = True
    if parser:
        parser = argparse.ArgumentParser("Generates conditional samples of system via Algabraic Decision Diagrams.")
        parser.add_argument("fname", help="Model exported as drdd file by storm", type=str)
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
        filename = "dtmcs/brp/brp_N_16_MAX_4.drdd"
        path_n = 64
        repeats = 100
        tlabel = 'target'
        store = False
    print(f'Running parameters: fname={filename}, n={path_n}, repeats={repeats}, label={tlabel}, store={store}')
    
    manager = _agd.ADD()
    manager.configure(max_growth=1.5)
    context = lambda: None # (required to assign attributes)
    context.manager = manager # type: ignore
    
    parse_time = perf_counter_ns()
    model = load_adds_from_drdd(context.manager, filename, # type: ignore
                                load_targets=['initial', 'transitions', f'label {tlabel}'])
    print(f'Finished parsing input: {ms_str_from(parse_time)}.')
    init = model['initial']
    target = model[f'label {tlabel}']
    assert len(target) > 0, "Target states missing"
    transitions = model['transitions']
    
    print(f"Number of variables per state: {len(transitions.support)//2}")
    print(f"Size of ADD: {transitions.dag_size} nodes")

    if store:
        raise NotImplementedError("ADD storage is not yet supported")
        # dirname = filename.replace('.drdd', '/')
        # gs, ts = load_and_store(dirname, transitions, path_n)
    else:
        precomp_time = perf_counter_ns()
        gs, ts = compute_power_graphs(context, transitions, path_n)
        print(f'Finished precomputing functions: {ms_str_from(precomp_time)}.')
        
    res = generate_many_traces(context, ts, path_n,
                init, target, save_traces=False,
                repeats=repeats)
    
