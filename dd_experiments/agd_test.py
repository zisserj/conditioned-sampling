import dd.cudd_add as _agd

def example_agd():
    agd = _agd.ADD()
    agd.declare('x', 'y', 'z')
    one = agd.constant(1)
    two = agd.constant(2)
    x = agd.var('x')
    y = agd.var('y')
    u = agd.apply('+', one, x)
    v = agd.apply('*', y, two)
    r = agd.apply('-', u, v)
    print(r)
    values = dict(
        x=True,
        y=True)
    p = agd.let(values, r)
    print(p)
    agd.dump('rooted.pdf', roots=[p])
    assert p == agd.constant(0)
    
    
def test_to_expr():
    agd = _agd.ADD()
    agd.declare('x', 'y')
    # x
    u = agd.var('x')
    expr = agd.to_expr(u)
    assert expr == 'x', expr
    # y
    u = agd.var('y')
    expr = agd.to_expr(u)
    assert expr == 'y', expr
    # x /\ ~ y
    u = agd.add_expr(r'x /\ ~ y')
    expr = agd.to_expr(u)
    expr_ = 'ite(x, ite(y, 0.0, 1.0), 0.0)'
    assert expr == expr_, expr
    print('completed testing `to_expr`')


def variable_substitution():
    # instantiate a shared BDD manager
    agd = _agd.ADD()
    agd.declare('x', 'y', 'u', 'v')
    # create the BDD for the disjunction of x and y
    u = agd.add_expr(r'x \/ y')
    # Substitution of x' for x and y' for y.
    # In TLA+ we can write this as:
    #
    # LET
    #     x == u
    #     y == v
    # IN
    #     x \/ y
    rename = dict(x='u', y='v')
    v = agd.let(rename, u)
    # show the result
    s = agd.to_expr(v)
    # agd.dump('rooted.pdf', roots=[v])
    print(s)

    # another way to confirm that the result is as expected
    v_ = agd.add_expr(r'u \/ v')
    assert v == v_

if __name__ == '__main__':
    example_agd()
    test_to_expr()
    variable_substitution()