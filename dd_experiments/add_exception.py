import dd.cudd_add as _agd

manager = _agd.ADD()

var_names = ['x', 'y', 'z'] # no exception if 'x' is removed
manager.declare(*var_names)
expr = manager.add_expr('y')

map_exists = {'y': 'y', 'z': 'z'} # no exception if either mapping is removed
w = manager.let(map_exists, expr)
