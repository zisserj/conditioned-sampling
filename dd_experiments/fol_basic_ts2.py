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
    ver_num = 4
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
        x=(0,ver_num-1),
        y=(0,ver_num-1),
        c0=(0,1),
        g0='bool')
    g = ctx.add_expr(exprs_union)
    val_cond = ctx.add_expr(r'(c0=1 <=> g0)')
    ret = ctx.replace_with_bdd(val_cond, {'g0': g})
    #print(ctx.to_expr(ret))
    return ret, ver_num

g0, ver_num = make_basic_ts_bdd()
#bdd.dump('dd_experiments/g0.pdf', [g0])

def mul_g0(g):
    ctx.declare(
        z=(0,ver_num-1),
        c0_=(0,1),
        g0_='bool',
        hc0=(0,1)) # half-count

    rename = dict(c0='c0_', x='y', y='z')

    g_ = ctx.let(rename, g)

    mult_cond = ctx.add_expr(r'g0 /\ g0_ /\ (hc0 = c0 * c0_)')
    mult_cond = ctx.replace_with_bdd(mult_cond, {'g0': g, 'g0_': g_})
    return ctx.exist({'c0', 'c0_'}, mult_cond)

def mul_gi(context, gi, i, g_c_max=1):
    vrs = {'z': (0,ver_num-1),
           f'c{i}_': (0, g_c_max),
           f'g{i}_': 'bool',
           f'hc{i}': (0, g_c_max**2)} #halfcount
    context.declare(**vrs)
    print(vrs)
    rename = {'x': 'y',
              'y': 'z',
           f'c{i}': f'c{i}_'}
    gi_ = context.let(rename, gi)

    mult_cond = context.add_expr(f'g{i} /\\ g{i}_ /\\ (hc{i} = c{i} * c{i}_)')
    mult_cond = context.replace_with_bdd(mult_cond, {f'g{i}': gi, f'g{i}_': gi_})
    return context.exist({f'c{i}', f'c{i}_'}, mult_cond)


# print(ctx.support(h_tag))
# print(bdd.support(h_tag))

def sum_to_g(t, g_c_max, context):
    #hs = halfstep
    hs_vars = [f't_{i}' for i in range(ver_num)]
    hs_counts = [f'val_{i}' for i in range(ver_num)]
    context.declare(**{h: 'bool' for h in hs_vars})
    context.declare(**{val: (0, g_c_max**2) for val in hs_counts})
    context.declare(ch_sum=(0, (g_c_max**2)*ver_num))

    # find half-count var name
    c_name = list(filter(lambda s: s.startswith('hc'), context.vars.keys()))[0]
    hs_y = [context.let({c_name:hc}, t) for hc in hs_counts]
    hs_bdds = [context.let(dict(y=i), ybdd) for i, ybdd in enumerate(hs_y)]
    hs_ti_expr = '/\\'.join(hs_vars)
    hs_sum = 'ch_sum = {}'.format('+'.join(hs_counts))
    hs_eq_vars = context.add_expr(f'{hs_ti_expr}  /\\ ({hs_sum})')
    hs_eq_full = context.replace_with_bdd(hs_eq_vars, {ch: cbdd for (ch, cbdd) in zip(hs_vars, hs_bdds)})
    g_ = context.exist(set(hs_counts), hs_eq_full)
    return g_


def make_next_ctx(context, i, gi_source, max_c_val):
    new_ctx = _fol.Context()
    new_ctx.declare(x=(0,ver_num-1),
        y=(0,ver_num-1),
        z=(0, ver_num-1))
    rename_vars = {f'c{i}': (0, max_c_val), f'g{i}': 'bool'}
    new_ctx.declare(**rename_vars)
    context.declare(**rename_vars)

    rename_gi = context.let({'ch_sum': f'c{i}', 'z':'y'}, gi_source)
    gi = context.copy(rename_gi, new_ctx)
    return (new_ctx, gi)


# (1) find highest bit that results in true valuation of BDD
# or (2) go over all true valuations of BDD to find largest c val
# this implements (1)


def simplfy_gs(context, g):
  i = context.i
  c_name = f'c{i}'
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



t0 = mul_g0(g0)
bdd.dump('t0.png', [t0])


ctxs = [ctx]
gs = [g0]
ts = []

gc_max = 1
ctx_k = ctx
g_k = g0
last_t = time_ns()
for k in range(0, 7):
    # t = g x g
    t_k = mul_gi(ctx_k, g_k, k, gc_max)
    ts.append(t_k)

    # sum t over y
    pre_g_k = sum_to_g(t_k, gc_max, ctx_k)

    gc_max = (gc_max**2) * ver_num

    # rename vars, make new ctx
    ctx_k, g_k = make_next_ctx(ctx_k, k+1, pre_g_k, gc_max)
    ctx_k.i = k+1 # type: ignore
    ctxs.append(ctx_k)
    
    print(f'Generated {k+1}: {time_ns()-last_t}')
    last_t = time_ns()
    # simplify varcount of g
    g_k = simplfy_gs(ctx_k, g_k)
    
    print(f'Simplified {k+1}: {time_ns()-last_t}')
    last_t = time_ns()
    
    gc_max = ctx_k.vars[f'c{k+1}']['dom'][1]
    gs.append(g_k)