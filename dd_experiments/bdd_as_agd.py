"""How to copy a BDD from one manager to another."""
import dd.cudd as _bdd
import dd.cudd_add as _agd


def transfer():
    """Copy a BDD from one manager to another."""
    # create two BDD managers
    source = _bdd.BDD()
    target = _agd.ADD()
    # declare the variables in both BDD managers
    vrs = ['a', 'b']
    source.declare(*vrs)
    target.declare(*vrs)
    # create a BDD with root `u`
    u = source.add_expr(r'a /\ b')
    print(source.vars)
    # copy the BDD `u` to the BDD manager `target`
    u_ = source.copy(u, target)



def copy_variable_order():
    """As in `transfer`, and copy variable order too."""
    source = _bdd.BDD()
    target = _agd.ADD()
    # declare variables in the source BDD manager
    source.declare('a', 'b')
    # create a BDD with root `u`
    u = source.add_expr(r'a /\ b')
    # copy the variables, and the variable order
    target.declare(*source.vars)
    target.reorder(source.var_levels)
    # copy the BDD `u` to the BDD manager `target`
    u_ = source.copy(u, target)

def dump_bdd_as_json(filename):
    """Write a BDD to a JSON file."""
    bdd = _bdd.BDD()
    bdd.declare('x', 'y', 'z')
    u = bdd.add_expr(r'(x /\ y) \/ ~ z')
    roots = dict(u=u)
    bdd.dump(filename, roots)
    print(f'Dumped BDD: {u}')

if __name__ == '__main__':
    transfer()
    copy_variable_order()