import dd.cudd_add as _agd # type: ignore
import itertools

print_refs = lambda f: print(f'attr: {f._ref}, cudd: {f.ref}')

def sus_func(o=None):
    manager = _agd.ADD()
    var_names = ['x', 'y', 'z'] # no exception if 'x' is removed (Cudd_addIthVar)
    manager.declare(*var_names)
    expr = manager.constant(2)

    if o:
        mymap = {o[0]:o[1], o[2]:o[3]}
    else:
        mymap = {'x': 'y', 'y': 'z'} # no exception if either is removed, bc its in multi_compose
    w1 = manager.let(mymap, expr)
    #print_refs(manager.constant(4))

opts = itertools.product('xyz', repeat=4)
# for o in opts:
#     if (o[0] != o[2]):
#         print(''.join(o))
#         sus_func(o)
#         print('-')
sus_func()