import dd.cudd_add as _agd
import math
from time import time_ns

from add_from_drdd import load_adds_from_drdd


def make_sample_add(manager):
    two = manager.constant(2)

    ver_num = 4
    req_vars = math.ceil(math.log2(ver_num))
    x_var_names = [f'x{i}' for i in range(req_vars)]
    manager.declare(*x_var_names)
    y_var_names = [f'y{i}' for i in range(req_vars)]
    manager.declare(*y_var_names)
    z_var_names = [f'z{i}' for i in range(req_vars)]
    manager.declare(*z_var_names)
    x_vars = [manager.var(xi) for xi in x_var_names]
    y_vars = [manager.var(yi) for yi in y_var_names]
    opr_x = manager.apply('+', *x_vars) # doesnt scale for more vars but later problem
    opr_y = manager.apply('*', *y_vars)
    #opr_y = manager.apply('+', opr_y, two)
    opr_xy = manager.apply('\\/', opr_x, opr_y)
    return opr_xy

def make_min_sample(manager):
    var_names = ['x1', 'x2', 'y1', 'y2']
    manager.declare(*var_names)
    u = manager.add_expr('x1 | x2 | y1 | y2')
    return u

def process_count_structure(manager, g0, maxN):
    # g0 must have x1-xn, y1-yn vars

    # find/define variables req
    x_var_names = []
    y_var_names = []
    for v in g0.support:
        if v.startswith('x'):
            x_var_names.append(v)
        elif v.startswith('y'):
            y_var_names.append(v)
        else:
            raise Exception("Unrecognised variable in g0")
    assert len(x_var_names) == len(y_var_names)
    x_var_names = sorted(x_var_names)
    y_var_names = sorted(y_var_names)
    z_var_names = [v.replace('x', 'z') for v in x_var_names]
    manager.declare(*z_var_names)

    map_mul = dict(zip(x_var_names, y_var_names)) | dict(zip(y_var_names, z_var_names))
    map_exists = dict(zip(z_var_names, y_var_names))

    gs = [g0]
    ts = []

    g_k = g0
    last_t = time_ns()
    for k in range(0, int(math.log2(maxN))):
        # t = g x g
        g_k_ = manager.let(map_mul, g_k)
        t_k = manager.apply('*', g_k, g_k_) # cuddGarbageCollect?
        ts.append(t_k)
        # g = Ey in t
        g_k_pre = manager.exist(y_var_names, t_k)
        g_k = manager.let(map_exists, g_k_pre)
        gs.append(g_k)
        
        print(f'Generated {k+1}: {(time_ns()-last_t)*1e-9}')
        last_t = time_ns()

    return gs, ts


manager = _agd.ADD()
#g0 = make_min_sample(manager)
filename = "dd_experiments/die.drdd"
#filename = "dtmcs/brp/brp_16_2.drdd"
adds = load_adds_from_drdd(manager, filename,
                                rename_vars=True)
    
g0 = adds['transitions']
gs, ts = process_count_structure(manager, g0, 8)
print('done')