import dd.cudd as _bdd
import dd.cudd_add as _agd

def dump_bdd_as_json(filename):
    """Write a BDD to a JSON file."""
    bdd = _bdd.BDD()
    bdd.declare('x', 'y', 'z')
    u = bdd.add_expr(r'(x /\ y) \/ ~ z')
    roots = dict(u=u)
    bdd.dump(filename, roots)
    print(f'Dumped BDD: {u}')


def load_add_from_json(filename):
    """Load a BDD from a JSON file."""
    agd = _agd.ADD()
    roots = agd.load(filename)
    print(f'Loaded AGD: {roots}')

if __name__ == '__main__':
    filename = 'storage.json'
    dump_bdd_as_json(filename)
    load_add_from_json(filename)