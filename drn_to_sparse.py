import re
import numpy as np
import scipy.sparse as sp

# storm --prism crowds.pm --constants "TotalRuns=3,CrowdSize=10" --buildfull --prismcompat --engine sparse --exportbuild test_base.drn
'''
storm --prism brp.pm --constants N=16,MAX=2 --buildfull --prismcompat --engine sparse --exportbuild brp_16_2.drn
storm --prism brp.pm --constants N=16,MAX=2 --buildfull --prismcompat --engine dd --exportbuild brp_16_2.drdd
'''


def read_drn(filename, target_label='target'):
    with open(filename) as f:
        content = f.read()

    # currently does not capture rewards
    pat_with_sqr = r"state (\d+)(?: \[[\d., ]+\])? ?([\w ]*)\n\taction 0(?: \[[\d., ]+\])?\n([\s\d:.]*)"

    rows_match = re.finditer(pat_with_sqr, content)
    num_states = int(re.findall(r"@nr_states\n(\d+)",content)[0])
    
    t_mat = sp.dok_array((num_states, num_states))
    initial_states = []
    target_states = []
    for match in rows_match:
        num_states += 1
        idx = int(match.group(1).strip())
        labels = match.group(2).strip().split()
        if "init" in labels:
            initial_states.append(idx)
        if target_label in labels:
            target_states.append(idx)
        body = match.group(3).strip()
        ts_strs = re.finditer(r"(\d+) : ([\d.]+)", body)
        for ts_match in ts_strs:
            t = int(ts_match.group(1))
            p = float(ts_match.group(2))
            t_mat[idx, t] = p
            #print(f'{idx} - {t}:{p} ({labels})')

    initial = np.array(initial_states)
    target = np.array(target_states)
    transitions = t_mat.tocsr()
    assert transitions.nnz > 0, "Input drn was not processed correctly"

    return {'init': initial, 'target': target, 'trans': transitions}

if __name__ == '__main__':
    filename = "dtmcs/brp/brp_64_2.drn"
    res = read_drn(filename)
    print(res['initial'])
    print(res['target'])
    # sp.save_npz("dice_mat.npz", transitions)
