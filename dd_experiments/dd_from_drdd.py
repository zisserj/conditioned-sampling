from fileinput import filename
import dd.cudd_add as _agd
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


def build_add(agd, vars, nodes_iter):
    agd.configure(reordering=False)
    agd.declare(*vars)
    temp_cache = [0]
    for s in nodes_iter:
        if s.group(1) == 'leaf':
            val = float(s.group(3))
            leaf = agd.constant(val)
            temp_cache.append(leaf)
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
            agd.declare(var)
            v = agd.var(var)
            u = agd.ite(v, temp_cache[high], temp_cache[low])
            temp_cache.append(u)
    res = temp_cache[-1]
    del temp_cache
    return res

def load_adds_from_drdd(agd, filename):
    with open(filename, "r") as file:
        content = file.read()

    adds_str = re.finditer(r"%([\S ]+)\n\[\n([\S\s]*?)\n\],\[(\d+),\]", content)
    res = {}
    for match in adds_str:
        name = match.group(1).strip()
        body = match.group(2).strip()
        nodes_iter = re.finditer(r"(leaf)\((\d+),\d+,\"([\d.]+)\"\)|(node)\((\d+),(\d+),(\d+),(~?\d+)\)", body)
    
        # need to have all variables declared in advance for some reason
        vars = re.findall(r"node\(\d+,(\d+),\d+,~?\d+\)", body)
        vars = sorted(list(set(vars)))
        vars = [f'v{i}' for i in vars]
        
        size = int(match.group(3).strip())
        res[name] = build_add(agd, vars, nodes_iter)
    return res
    

if __name__ == '__main__':
    agd = _agd.ADD()

    filename = "/home/jules/storm_sampler/storm-project-starter-cpp/symbolic_model.drdd"
    # filename = 'dd_experiments/test.drdd'
    adds = load_adds_from_drdd(agd, filename)
    for name, add in adds.items():
        print(f'{name} has {add.dag_size} nodes')
        print(f'support = {add.support}')