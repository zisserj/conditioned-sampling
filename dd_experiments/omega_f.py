import omega.symbolic.fol as _fol
#import dd.cudd as _bdd
_bdd = _fol._bdd

ctx = _fol.Context()
ctx.declare(
    x='bool',
    y=(0, 5))
u: _bdd.Function = ctx.add_expr(
    r' x /\ (y = 5) ')


support = ctx.support(u)
assert support == {'x', 'y'}
support_bits = ctx.bdd.support(u)
assert support_bits == {'x', 'y_0', 'y_1', 'y_2'}

ctx.bdd.dump('dd_experiments/test.pdf', [u])