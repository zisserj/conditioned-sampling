from fileinput import filename
import dd.cudd as _bdd
import omega.symbolic.fol as _fol
import re



''' https://github.com/trolando/sylvan/blob/master/src/sylvan_mtbdd.h#L853
 * Write <count> decision diagrams given in <dds> in ASCII form to <file>.
 * Also supports custom leaves using the leaf_to_str callback.
 *
 * The text format writes in the same order as the binary format, except...
 * [
 *   node(id, var, low, high), -- for a normal node (no complement on high)
 *   node(id, var, low, ~high), -- for a normal node (complement on high)
 *   leaf(id, type, "value"), -- for a leaf (with value between "")
 * ],[dd1, dd2, dd3, ...,] -- and each the stored decision diagram.
'''


def build_bdd(ctx, vars, nodes_iter, prob=True, denominator=1000):
    # probability denominator is implicit in the bdd 
    ctx.bdd.configure(reordering=False)
    ctx.declare(**{v: 'bool' for v in vars})
    if prob:
        ctx.denominator = denominator
        ctx.declare(p0=(0,denominator))
    temp_cache = [None]
    for s in nodes_iter:
        if s.group(1) == 'leaf':
            val = float(s.group(3))
            
             #v1< ((p * p_)/(dnm))) < v2
            
            
            if prob:
                rounded_val = int(val * denominator)
                node = ctx.add_expr(f'p0={rounded_val}')
            else:
                node = ctx.true if val == 1 else ctx.false
            temp_cache.append(node) # type: ignore
        elif s.group(4) == 'node':
            var ='v' + s.group(6)
            low = int(s.group(7))
            c = False
            if '~' in s.group(8):
                high = int(s.group(8)[1:])
                c = True
                raise ValueError("Complemented high nodes are not supported (yet)")
            else:
                high = int(s.group(8))
            v = ctx.to_bdd(var)
            u = ctx.add_expr(f"IF {v} THEN {temp_cache[high]} ELSE {temp_cache[low]}")
            temp_cache.append(u) # type: ignore
    res = temp_cache[-1]
    del temp_cache
    ctx.bdd.configure(reordering=True)
    return res

def rename_vars_xy(ctx, bdd_dict, vars):
    # assuming storm export, where if v(x)=d then v(x+1)=d'
    map = {}
    for i, v in enumerate(vars):
        new_var = 'x' if i % 2 == 0 else 'y'
        map[v] = f'{new_var}_{i // 2}'
    vert_domain = (0, 2**(len(vars)//2))
    ctx.declare(**{v: vert_domain for v in ['x', 'y', 'z']})
    
    for g_name in bdd_dict.keys():
        g = bdd_dict[g_name]
        bdd_dict[g_name] = ctx.bdd.let(map, g)
    ctx.vars = {k:v for k, v in ctx.vars.items() if k in ['x', 'y', 'z', 'p0']}

def load_bdds_from_drdd(ctx, filename,
                        rename_vars=True, load_targets=['transitions'],
                        denominator=1000):
    with open(filename, "r") as file:
        content = file.read()

    vars = []
    adds_str = re.finditer(r"%([\S ]+)\n\[\n([\S\s]*?)\n\],\[(\d+),\]", content)
    res = {}
    for match in adds_str:
        name = match.group(1).strip()
        body = match.group(2).strip()
        nodes_iter = re.finditer(r"(leaf)\((\d+),\d+,\"([\d.]+)\"\)|(node)\((\d+),(\d+),(\d+),(~?\d+)\)", body)
    
        # need to have all variables declared in advance or ADD levels get messed up
        if name == 'transitions':
            vars = re.findall(r"node\(\d+,(\d+),\d+,~?\d+\)", body)
            vars = sorted(list(set(vars)))
            vars = [f'v{i}' for i in vars]

        size = int(match.group(3).strip())
        if name in load_targets:
            if name == "transitions":
                res[name] = build_bdd(ctx, vars, nodes_iter, True, denominator)
            else:
                res[name] = build_bdd(ctx, vars, nodes_iter, False)
    
    if rename_vars:
        rename_vars_xy(ctx, res, vars)
    return res


if __name__ == '__main__':
    bdd_manager = _bdd.BDD()
    ctx = _fol.Context()
    ctx.bdd = bdd_manager # use cudd instead of python impl


    #filename = "/home/jules/storm_sampler/storm-project-starter-cpp/symbolic_model.drdd"
    filename = "dd_experiments/die.drdd"
    targets = ['transitions', 'initial', 'label target', 'label one']
    bdds = load_bdds_from_drdd(ctx, filename, load_targets=targets)

    # change everything but prob to use 0/1
    for name, bdd in bdds.items():
        print(f"{name}'s support is {ctx.support(bdd)}")


    