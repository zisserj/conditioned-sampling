import re
import numpy as np


fname = 'dtmcs/dice/die.drn.out'
with open(fname) as f:
    content = f.read()
labels_content, traces_content = content.split('---')

matches = re.finditer(r'([\w]+): \[([\d+,\s]+)\]', labels_content)

labels_dict = {}
for m in matches:
    labels_dict[m.group(1)] = np.fromstring(m.group(2), sep=' ', dtype=int)


n_heads = 0
n_tails = 0
total_gr_heads = 0

traces = traces_content.split('\n')[1:]
for t in traces:
    t_arr = np.fromstring(t, sep=',', dtype=int)
    # state_visited = np.isin(t_arr, labels_dict['all_visited'])
    check_heads = np.isin(t_arr, labels_dict['flip_heads'])
    #finished_idx = np.searchsorted(state_visited,True)
    n_heads += check_heads.sum()
    
    check_tails = np.isin(t_arr, labels_dict['flip_tails'])
    #finished_idx = np.searchsorted(state_visited,True)
    n_tails += check_tails.sum()
    # if t_arr[finished_idx] in labels_dict['on_hand6']:
    #     count_on6 += 1
    #     total += finished_idx
    if check_heads.sum() >= check_tails.sum():
        total_gr_heads += 1
print('avg throws: ', (n_heads+n_tails)/len(traces))
print('avg heads: ', n_heads/len(traces))
print('avg tails: ', n_tails/len(traces))
print('P[n_h>n_t]: ', total_gr_heads/len(traces))