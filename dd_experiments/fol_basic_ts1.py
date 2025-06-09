import omega.symbolic.fol as _fol
#import dd.cudd as _bdd
_bdd = _fol._bdd

ctx = _fol.Context()
bdd = ctx.bdd

ctx.declare(
    c=(0,1),
    x=(1,3),
    y=(1,3))

base = r'((c={}) /\ (x={}) /\ (y={}))'
trans = '''1, 1, 1
1, 1, 2
1, 2, 3
1, 3, 3
0, 1, 3
0, 2, 1
0, 2, 2
0, 3, 1
0, 3, 2
'''.replace(' ', '')

ts = [l.split(',') for l in trans.split()]
exprs = [base.format(*vars) for vars in ts]
union = r'\/'.join(exprs)
# u: _bdd.Function = ctx.add_expr(
#     r'(c=1) /\ (x=1) /\ (y=2)')
u: _bdd.Function = ctx.add_expr(union)

support = ctx.support(u)
#assert support == {'x', 'y'}
support_bits = bdd.support(u)
#assert support_bits == {'x', 'y_0', 'y_1', 'y_2'}

bdd.dump('dd_experiments/test.pdf', [u])

ctx.declare(
    c2=(0,1),
    z=(1,3)
)

rename = dict(c='c2', x='y', y='z')
u2 = ctx.let(rename, u)
bdd.dump('dd_experiments/test1.pdf', [u2])


ctx.declare(
    next_c=(0,3)
)

print(ctx.to_expr(u,))
#print(ctx.to_expr(u2))