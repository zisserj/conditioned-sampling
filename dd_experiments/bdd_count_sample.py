import omega.symbolic.fol as _fol
import omega.symbolic.functions as _fun
from time import time_ns

_bdd = _fol._bdd


# everything using the same global context for now
ctx = _fol.Context()
bdd = ctx.bdd
#bdd.configure(reordering=True) # already the default in cudd


def make_basic_ts_bdd():
    r"""
    recieves fol context, defines
    g0(x, y, c0) <=> g(x, y) = c0
    returns:
    g0 bdd, var_num vertices in ts
    """
    base = r'((x={}) /\ (y={}))'
    vert_num = 4
    trans = '''0, 1
    0, 3
    1, 1
    1, 2
    2, 2
    3, 1'''.replace(' ', '')
    ts = [l.split(',') for l in trans.split()]
    exprs = [base.format(*vars) for vars in ts]
    exprs_union = r'\/'.join(exprs)

    ctx.declare(
        x=(0,vert_num-1),
        y=(0,vert_num-1),
        z=(0,vert_num-1),
        c0=(0,1),
        g='bool',
        g_='bool')
    g = ctx.add_expr(exprs_union)
    val_cond = ctx.add_expr(r'(c0=1 <=> g)')
    ret = ctx.replace_with_bdd(val_cond, {'g': g})
    #print(ctx.to_expr(ret))
    return ret, vert_num

g0, vert_num = make_basic_ts_bdd()
#bdd.dump('dd_experiments/g0.pdf', [g0])


def mul_gi(context, gi, i, g_c_max=1):
    vrs = {f'c{i}_': (0, g_c_max),
           f'mid{i}': (0, g_c_max**2),
           'g':'bool',
           'g_':'bool'} #halfcount
    context.declare(**vrs)
    rename = {'x': 'y',
              'y': 'z',
           f'c{i}': f'c{i}_'}
    gi_ = context.let(rename, gi)

    mult_cond = context.add_expr(f'g /\\ g_ /\\ (mid{i} = c{i} * c{i}_)')
    mult_cond = context.replace_with_bdd(mult_cond, {'g': gi, 'g_': gi_})
    exist = context.exist({f'c{i}', f'c{i}_'}, mult_cond) 
    context.vars.pop(f'c{i}', None) # remove old c var
    return exist


# print(ctx.support(h_tag))
# print(bdd.support(h_tag))

def sum_to_g(context, t, i, g_c_max):
    #hs = halfstep
    
    # identical every iteration
    ts_vars = [f't_{j}' for j in range(vert_num)]
    context.declare(**{h: 'bool' for h in ts_vars})
    hs_ti_expr = '&'.join(ts_vars)
    
    hs_counts = [f'val_{j}' for j in range(vert_num)]
    context.declare(**{val: (0, g_c_max**2) for val in hs_counts})
    
    context.declare(**{f'sum{i}':(0, (g_c_max**2)*vert_num)})

    hs_y = [context.let({f'mid{i}':hc}, t) for hc in hs_counts]
    hs_bdds = [context.let(dict(y=i), ybdd) for i, ybdd in enumerate(hs_y)]
    
    hs_sum = 'sum{} = {}'.format(i, '+'.join(hs_counts))
    hs_eq_vars = context.add_expr(f'{hs_ti_expr}  /\\ ({hs_sum})')
    hs_eq_full = context.replace_with_bdd(hs_eq_vars, {ch: cbdd for (ch, cbdd) in zip(ts_vars, hs_bdds)})
    g_ = context.exist(set(hs_counts), hs_eq_full)
    
    # remove placeholder variables
    [context.vars.pop(val, None) for val in hs_counts]
    return g_


def make_next_iter(context, next_i, gi, max_c_val):
    rename_vars = {f'c{next_i}': (0, max_c_val)}
    context.declare(**rename_vars)
    
    
    new_ctx = _fol.Context()
    new_ctx.declare(x=(0,vert_num-1),
                    y=(0,vert_num-1),
                    z=(0, vert_num-1))
    new_ctx.declare(**rename_vars)
    
    new_g = context.let({f'sum{next_i-1}': f'c{next_i}', 'z':'y'}, gi)
    gi = context.copy(new_g, new_ctx)
    return (new_ctx, gi)


# (1) find highest bit that results in true valuation of BDD
# or (2) go over all true valuations of BDD to find largest c val
# this implements (1)


def simplfy_gs(context, k, g):
  c_name = f'c{k}'
  bits = context.vars[c_name]['bitnames']

  bits_to_remove = []
  for c_bit in reversed(bits): # switch to binary search
      res = g.bdd.let({c_bit:True}, g)
      # possibly use faster check here?
      g.bdd.dump('c1_2_true.png', [res])
      if res.count() > 0:
          break
      bits_to_remove.append(c_bit)
  if len(bits_to_remove) == 0:
      return g
  g_trimmed = g.bdd.let({extra_bit: False for extra_bit in bits_to_remove}, g)
  new_width = len(bits) - len(bits_to_remove)
  context.vars[c_name]['bitnames'] = bits[:new_width]
  context.vars[c_name]['width'] = new_width
  context.vars[c_name]['dom'] = (0, 2**new_width -1)
  return g_trimmed



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