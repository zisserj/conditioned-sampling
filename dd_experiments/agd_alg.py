import dd.cudd_add as _agd
import math
from time import time_ns

import omega.symbolic.fol as _fol

manager = _agd.ADD()
manager.declare('w', 'm')
# create the BDD for the disjunction of x and y
zero = manager.constant(0)
one = manager.constant(1)
two = manager.constant(2)
w = manager.var('w')
m = manager.var('m')
u = manager.apply('+', one, w)
v = manager.apply('*', m, two)
r = manager.apply('-', u, v)

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
z_vars = [manager.var(zi) for zi in z_var_names]

opr_x = manager.apply('+', *x_vars) # doesnt scale for more vars but later problem
opr_y = manager.apply('*', *y_vars)
opr_y = manager.apply('+', opr_y, two)
opr_xy = manager.apply('\\/', opr_x, opr_y)
manager.dump('dd_experiments/opr_all.png', [opr_xy])



rename_mul = dict(zip(x_var_names, y_var_names)) | dict(zip(y_var_names, z_var_names))
rename_exists = dict(zip(y_var_names, z_var_names))

g0 = opr_xy
g0_ = manager.let(rename_mul, g0)
t0 = manager.apply('*', g0, g0_)
#manager.dump('dd_experiments/g0_mul.png', [t0])

g1 = manager.exist(y_var_names, t0)
#manager.dump('dd_experiments/g1.png', [g1])


gs = [g0]
ts = []

g_k = g0
last_t = time_ns()
for k in range(0, 7):
    # t = g x g
    g_k_ = manager.let(rename_mul, g_k)
    t_k = manager.apply('*', g_k, g_k_)

    ts.append(t_k)

    # rename vars
    g_k_pre = manager.exist(y_var_names, t_k)
    g_k = manager.let(rename_exists, g_k_pre)
    
    gs.append(g_k)
    
    print(f'Generated {k+1}: {time_ns()-last_t}')
    last_t = time_ns()

print(manager.to_expr(g_k))