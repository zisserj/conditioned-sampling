import dd.cudd as _bdd
import dd.cudd_add as _agd

def dump_bdd_as_json(filename):
    """Write a BDD to a JSON file."""
    bdd = _bdd.BDD()
    bdd.declare('x', 'y', 'z', 'a', 'b')
    u = bdd.add_expr(r'(x /\ y) \/ ~ z')
    v = bdd.add_expr(r'(a /\ b) \/ ~ (x \/ y)')
    w = u & v
    bdd.dump(filename, [w])
    bdd.dump(filename.replace('json','png'), [w])
    print(f'Dumped BDD: {u}')

def dump_add_as_json(filename):
    """Write a ADD to a JSON file."""
    agd = _agd.ADD()
    agd.declare('x', 'y', 'z')
    one = agd.constant(1)
    two = agd.constant(2)
    x = agd.var('x')
    y = agd.var('y')
    u = agd.apply('+', one, x)
    v = agd.apply('*', y, two)
    r = agd.apply('-', u, v)
    # values = dict(
    #     x=True,
    #     y=True)
    # p = agd.let(values, r)
    # print(p)
    agd.dump(filename.replace('json','png'), [r])
    print(f'Dumped ADD: {r}')



def load_bdd_from_json(filename):
    """Load a BDD from a JSON file."""
    bdd = _bdd.BDD()
    roots = bdd.load(filename)
    print(f'Loaded BDD: {roots}')

def load_add_from_json(filename):
    """Load a ADD from a JSON file."""
    agd = _agd.ADD()
    roots = agd.load(filename)
    print(f'Loaded AGD: {roots}')

if __name__ == '__main__':
    filename = 'dd_experiments/bdd_storage.json'
    dump_bdd_as_json(filename)
    load_bdd_from_json(filename)