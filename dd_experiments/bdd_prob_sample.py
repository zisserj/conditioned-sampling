import omega.symbolic.fol as _fol
import dd.cudd as _bdd
from bdd_from_drdd import load_bdds_from_drdd
from time import time_ns




bdd = _bdd.BDD()
ctx = _fol.Context()
ctx.bdd = bdd

INITIAL_frac=100

file_bdds = load_bdds_from_drdd(ctx, "dd_experiments/die.drdd",
                                load_targets=['initial', 'label target', 'transitions'],
                                denominator=INITIAL_frac)
g0 = file_bdds['transitions']

def mul_gi(context, gi, i):
    dnm = context.vars[f'p{i}']['dom'][1] # not a tight bound
    vrs = {f'p{i}_': (0, dnm),
           f'mul_p{i}': (0, dnm),
           'g':'bool',
           'g_':'bool'}
    context.declare(**vrs)
    rename = {'x': 'y',
              'y': 'z',
           f'p{i}': f'p{i}_'}
    gi_ = context.let(rename, gi)
    lower_b = context.add_expr(f'(mul_p{i} * {dnm} - {dnm//2})< (p{i} * p{i}_)') # <= (r_p*{dnm} + {dnm//2})')
    upper_b = context.add_expr(f'(p{i} * p{i}_) <= (mul_p{i}*{dnm} + {dnm//2})')
    mult_g = gi & gi_ & lower_b & upper_b
    exist = context.exist({f'p{i}', f'p{i}_'}, mult_g)
    return exist


t0 = mul_gi(ctx, g0, 0)
gen_diff_mids = ctx.exist(['x', 'y', 'z'], t0)
for i in ctx.pick_iter(gen_diff_mids, ['mul_p0']):
    print(i)
# print(ctx.support(h_tag))
# print(bdd.support(h_tag))

def sum_to_g(context, t, i):
    #hs = halfstep
    
    vert_num = context.vars['x']['dom'][1]
    # identical every iteration
    ts_vars = [f't_{j}' for j in range(vert_num)]
    context.declare(**{h: 'bool' for h in ts_vars})
    hs_ti_expr = context.add_expr('&'.join(ts_vars))
    
    dnm = context.denominator
    hs_counts = [f'val_{j}' for j in range(vert_num)]
    context.declare(**{val: (0, dnm) for val in hs_counts})
    context.declare(**{f'sum{i}':(0, dnm*vert_num)})
    context.declare(**{f'psum{j}':(0, dnm*(j+1)) for j in range(vert_num)})

    # create t(x, 0, z, v0) & t(x, 1, z, v1) & ...
    hs_y = [context.let({f'mul_p{i}':hc}, t) for hc in hs_counts]
    hs_bdds = [context.let(dict(y=i), ybdd) for i, ybdd in enumerate(hs_y)]
    yi_transitions = context.replace_with_bdd(hs_ti_expr, {ch: cbdd for (ch, cbdd) in zip(ts_vars, hs_bdds)})
    
    
    context.add_expr('psum0 = val_0')
    cur_sum = yi_transitions
    for j in range(0, vert_num):
        if j == 0:
            cur_sum = cur_sum & context.add_expr('psum0 = val_0')
            cur_sum = context.exist({'val_0'} ,cur_sum)
        else:
            j_expr = context.add_expr(f'psum{j} = psum{j-1} + val_{j}')
            cur_sum = j_expr & cur_sum
            cur_sum = context.exist({f'psum{j-1}', f'val_{j}'} ,cur_sum) # don't care how current sum was reached            
    g_ = context.replace(cur_sum, {f'psum{vert_num-1}': f'sum{i}'})
    return g_

g1 = sum_to_g(ctx, t0, 0)

def make_next_iter(context, next_i, gi, max_c_val):
    rename_vars = {f'c{next_i}': (0, max_c_val)}
    context.declare(**rename_vars)
    vert_num = context.vars['x']['dom'][1]
    
    new_ctx = _fol.Context()
    new_ctx.declare(x=(0,vert_num-1),
                    y=(0,vert_num-1),
                    z=(0, vert_num-1))
    new_ctx.declare(**rename_vars)
    
    new_g = context.let({f'sum{next_i-1}': f'c{next_i}', 'z':'y'}, gi)
    gi = context.copy(new_g, new_ctx)
    return (new_ctx, gi)


ctxs = [ctx]
gs = [g0]
ts = []

gc_max = 1
ctx_k = ctx
g_k = g0
last_t = time_ns()
for k in range(0, 8):
    # t = g x g
    t_k = mul_gi(ctx_k, g_k, k, gc_max)
    ts.append(t_k)

    # sum t over y
    pre_g_k = sum_to_g(ctx_k, t_k, k, gc_max)

    gc_max = (gc_max**2) * vert_num

    # rename vars for next iter
    ctx_k, g_k = make_next_iter(ctx_k, k+1, pre_g_k, gc_max)
    ctxs.append(ctx_k)
    
    print(f'Generated {k+1}: {(time_ns()-last_t)*1e-9}')
    last_t = time_ns()
    # simplify varcount of g
    g_k = simplfy_gs(ctx_k, k+1, g_k)
    
    print(f'Simplified {k+1}: {(time_ns()-last_t)*1e-9}')
    last_t = time_ns()
    
    gc_max = ctx_k.vars[f'c{k+1}']['dom'][1]
    gs.append(g_k)