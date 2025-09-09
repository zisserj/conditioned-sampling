from time import perf_counter_ns
import omega.symbolic.fol as _fol
import dd.cudd as _bdd # type: ignore
from drdd_to_bdd import load_bdds_from_drdd
import numpy as np
#from line_profiler import profile

np.set_printoptions(precision=2, suppress=True)
rng = np.random.default_rng()
import argparse

ms_str_from = lambda start_ns: f'{(perf_counter_ns()-start_ns)*1e-6:05.6f}ms'
ms_str_any = lambda ns: f'{ns*1e-6:.6f}ms'


def compute_mid_step(ctx, gi, i):
    dnm = ctx.vars[f'p{i}']['dom'][1] # not a tight bound
    vrs = {f'p{i}_': (0, dnm),
           f'mul_p{i}': (0, dnm),
           'g':'bool',
           'g_':'bool'}
    ctx.declare(**vrs)
    rename = {'x': 'y',
              'y': 'z',
           f'p{i}': f'p{i}_'}
    gi_ = ctx.let(rename, gi)
    lower_b = ctx.add_expr(f'(mul_p{i} * {dnm} - {dnm//2})< (p{i} * p{i}_)') # <= (r_p*{dnm} + {dnm//2})')
    upper_b = ctx.add_expr(f'(p{i} * p{i}_) <= (mul_p{i}*{dnm} + {dnm//2})')
    mult_g = gi & gi_ & lower_b & upper_b
    exist = ctx.exist({f'p{i}', f'p{i}_'}, mult_g)
    return exist


def sum_to_g(ctx, t, i):
    # currently vert_num is upper bound from number of bits
    vert_num = ctx.vars['x']['dom'][1]
    
    dnm = ctx.denominator
    hs_counts = [f'val_{j}' for j in range(vert_num)]
    ctx.declare(**{val: (0, dnm) for val in hs_counts})
    ctx.declare(**{f'sum{i}':(0, dnm)})
    ctx.declare(**{f'psum{j}':(0, dnm) for j in range(vert_num)})

    # create t(x, 0, z, v0) & t(x, 1, z, v1) & ...
    yi_transitions = ctx.true
    for j in range(vert_num):
        bdd_j = ctx.let({f'mul_p{i}':f'val_{j}'}, t)
        bdd_j = ctx.let({'y':j}, bdd_j)
        yi_transitions &= bdd_j #& ctx.add_expr(f'val_{i} > 0'))
    
    # iterate on sums: psum(i) = psum(i-1) + val(i)
    cur_sum = yi_transitions
    for j in range(0, vert_num):
        if j == 0:
            cur_sum &= ctx.add_expr('psum0 = val_0')
            cur_sum = ctx.exist({'val_0'} ,cur_sum)
        else:
            j_expr = ctx.add_expr(f'psum{j} = psum{j-1} + val_{j}')
            cur_sum &= j_expr
            cur_sum = ctx.exist({f'psum{j-1}', f'val_{j}'} ,cur_sum) # don't care how current sum was reached            
    g_ = ctx.replace(cur_sum, {f'psum{vert_num-1}': f'sum{i}'})
    return g_

def sum_to_g_refined(ctx, t, i):
    dnm = ctx.denominator
    hs_counts = [f'val_{j}' for j in ctx.y_dom]
    ctx.declare(**{val: (0, dnm) for val in hs_counts})
    ctx.declare(**{f'sum{i}':(0, dnm)})
    ctx.declare(**{f'psum{j}':(0, dnm) for j in ctx.y_dom})

    # create t(x, 0, z, v0) & t(x, 1, z, v1) & ...
    yi_transitions = ctx.true
    for j in ctx.y_dom:
        bdd_j = ctx.let({f'mul_p{i}':f'val_{j}'}, t)
        bdd_j = ctx.let({'y':j}, bdd_j)
        yi_transitions &= bdd_j #& ctx.add_expr(f'val_{i} > 0'))
    
    # iterate on sums: psum(i) = psum(i-1) + val(i)
    cur_sum = yi_transitions
    last_j = -1
    for j in ctx.y_dom:
        if last_j == -1:
            cur_sum &= ctx.add_expr(f'psum{j} = val_{j}')
            cur_sum = ctx.exist({f'val_{j}'} ,cur_sum)
        else:
            j_expr = ctx.add_expr(f'psum{j} = psum{last_j} + val_{j}')
            cur_sum &= j_expr
            cur_sum = ctx.exist({f'psum{last_j}', f'val_{j}'} ,cur_sum) # don't care how current sum was reached            
        last_j = j
    g_ = ctx.replace(cur_sum, {f'psum{last_j}': f'sum{i}'})
    return g_

def make_next_iter_ctx(ctx, next_i, gi):
    rename_vars = {f'p{next_i}': ctx.vars[f'p{next_i-1}']['dom']}
    ctx.declare(**rename_vars)
    vert_num = ctx.vars['x']['dom'][1]
    
    new_ctx = _fol.Context()
    new_ctx.declare(x=(0,vert_num-1),
                    y=(0,vert_num-1),
                    z=(0, vert_num-1))
    new_ctx.declare(**rename_vars)
    new_ctx.denominator = ctx.denominator
    new_g = ctx.let({f'sum{next_i-1}': f'p{next_i}', 'z':'y'}, gi)
    gi = ctx.copy(new_g, new_ctx)
    return (new_ctx, gi)


def rename_iter_vars(ctx, next_i, gi):
    rename_vars = {f'p{next_i}': ctx.vars[f'p{next_i-1}']['dom']}
    ctx.declare(**rename_vars)
    new_g = ctx.let({f'sum{next_i-1}': f'p{next_i}', 'z':'y'}, gi)
    return new_g


def compute_power_graphs(ctx, trans, length):
    gs = [trans]
    ts = []
    g_k = trans

    last_t = perf_counter_ns()
    for i in range(0, int(np.log2(length))):
        # t = g x g
        t_k = compute_mid_step(ctx, g_k, i)
        ts.append(t_k)
        
        # no need to sum G if its the last iteration
        if i < int(np.log2(length))-1:
            # sum t over y
            pre_g_k = sum_to_g_refined(ctx, t_k, i)
            g_k = rename_iter_vars(ctx, i+1, pre_g_k)
            gs.append(g_k)
        
        # print(f'Finished iteration {i}: {(perf_counter_ns()-last_t)*1e-9}')
        last_t = perf_counter_ns()
    leave_vars = ['x','y','z']
    ctx.vars = {k:v for k,v in ctx.vars.items() if k in leave_vars or k.startswith('p') or k.startswith('mul_p')}
    return gs, ts

def weighted_sample(opts_iter, p_var):
    coords = []
    weights = []
    for o in opts_iter:
        coords.append([o['x'], o['y'], o['z']])
        weights.append(o[p_var])
    coords = np.array(coords)
    weights = np.array(weights,dtype=float)
    if len(coords) == 0:
        return "No matching traces"
    weights /= weights.sum()
    return rng.choice(coords, axis=0, p=weights)


def sample_bdd_conditioned(ctx, t, start, target, w):
    target_rename = ctx.let({'x': 'z'}, target)
    p_label = [k for k in ctx.support(t) if k.startswith('mul_p')][0]
    non_zero = ctx.add_expr(f'{p_label} > 0')
    rel_states = t & start & target_rename & non_zero
    opts = ctx.pick_iter(rel_states, ['x','y','z',p_label])
    res = weighted_sample(opts, p_label)
    if type(res) == str:
        return res
    w[0] = res[0]
    w[len(w)//2] = res[1]
    w[-1] = res[2]

def sample_bdd_seq(ctx, t, x_idx, z_idx, w):
    start_bdd = ctx.add_expr(f'x={w[x_idx]}')
    target_bdd = ctx.add_expr(f'z={w[z_idx]}')
    p_label = [k for k in ctx.support(t) if k.startswith('mul_p')][0]
    non_zero = ctx.add_expr(f'{p_label} > 0')
    rel_states = t & start_bdd & target_bdd & non_zero
    opts = ctx.pick_iter(rel_states, ['x','y','z',f'{p_label}'])
    res = weighted_sample(opts, p_label)
    w[(x_idx+z_idx)//2] = res[1]


def draw_sample(ctx, ts, length, init, target):
    w = np.full(length+1, -1, dtype=int)
    no_states = sample_bdd_conditioned(ctx, ts[-1], init, target, w)
    if no_states:
        return no_states
    for i in range(int(np.log2(length))-1, 0, -1):
        inc = np.power(2, i)
        for j in range(0, length, inc):
            sample_bdd_seq(ctx, ts[i-1], j, j + inc, w)
    return w
    
def state_to_og_vars(vars, total_bits, intval):
    bits = f'{intval:0{total_bits}b}'[::-1]
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

def print_gs(ctx, gs):
    no_zeros = [ctx.add_expr(f'p{i} > 0') for i in range(len(gs))]
    vars = [('s', 3), ('d', 3)]
    for g, nz in zip(gs, no_zeros):
        print('-----')
        g_ = g & nz
        asgns = list(ctx.pick_iter(g_))
        for a in asgns:
            x = state_to_og_vars(vars, 6, a['x'])
            y = state_to_og_vars(vars, 6, a['y'])
            print(f'x = {x}, y={y}')

if __name__ == "__main__":

    parser = True
    if parser:
        parser = argparse.ArgumentParser("Generates conditional samples of system via Boolean Decision Diagrams.")
        parser.add_argument("fname", help="Model exported as drdd file by storm", type=str)
        parser.add_argument("length", help="Generated trace length (currently only supports powers of 2)", type=int)
        parser.add_argument("precision", help="Number of bits used as denominator of all probabilies", type=int)
        parser.add_argument("-repeats", help="Number of traces to generate", type=int, default=1000)
        parser.add_argument("-tlabel", help="Name of target label matching desired final states",
                            type=str, default='target')
        parser.add_argument("-maxmem", help="Memory allocated to CUDD in GiB", type=int, default=1)
        parser.add_argument('--store', help="Store / try loading existing mats", action='store_true')
        args = parser.parse_args()
        filename = args.fname
        path_n = args.length
        precision = args.precision
        repeats = args.repeats
        tlabel = args.tlabel
        max_mem = args.maxmem
        store = args.store
    else:
        filename = "dtmcs/brp/brp_16_2.drdd"
        path_n = 4
        precision = 3
        repeats = 100
        tlabel = 'target'
        max_mem = 1
        store = False
    frac = 2**(precision) -1

    print(f'Running parameters: fname={filename}, n={path_n}, precision={precision},\
        repeats={repeats}, label={tlabel}, maxmem={max_mem}GiB, store={store}')
    
    bdd = _bdd.BDD()
    bdd.configure(max_memory = max_mem*(2**30))
    context = _fol.Context()
    context.bdd = bdd
    
    parse_time = perf_counter_ns()
    model = load_bdds_from_drdd(context, filename, # type: ignore
                                load_targets=['initial', 'transitions', f'label {tlabel}'],
                                denominator=frac)
    print(f'Finished parsing input: {ms_str_from(parse_time)}.')
    init = model['initial']
    target = model[f'label {tlabel}']
    assert len(target) > 0, "Target states missing"
    transitions = model['transitions']
    
    print(f"Number of variables per state: {context.vars['x']['width']}")
    print(f"Size of BDD: {transitions.dag_size} nodes")


    if store:
        raise NotImplementedError("BDD storage is not yet supported")
        # dirname = filename.replace('.drdd', '/')
        # gs, ts = load_and_store(dirname, transitions, path_n)
    else:
        precomp_time = perf_counter_ns()
        gs, ts = compute_power_graphs(context, transitions, path_n)
        print(f'Finished precomputing functions: {ms_str_from(precomp_time)}.')
        

    w = draw_sample(context, ts, path_n, init, target)
    res = generate_many_traces(context, ts, path_n,
                init, target, save_traces=True,
                repeats=10)
    

    # w = draw_sample(context, ts, path_n, init, target)
    # print(w)
    # vars = [('s', 3), ('d', 3)]
    
    # print([state_to_og_vars(vars, 6, wi) for wi in w])
