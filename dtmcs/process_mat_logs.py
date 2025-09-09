import re


pat = r"--- \w+\/([a-zA-Z]+)_?([\S_]+).drn - (\d+) ---\
[\s\S]*?\
[\w ]+input: ([\d.]+)ms.\
[\w ]+states: ([\d.]+)\
[\w ]+transitions: ([\d.]+)\
[\s\S]*?\
?(?:[\w ]+functions: ([\d.]+)ms.\
[\w ]+ mats: G\(\d+\), T\((\d+)\)|Found all required mats\.)\
[\w ]+ probability is ([\d.+-e]+)\
(?:[#\d\-\s]*Taken ([\d.]+)ms per sample|No matching traces)"



content = ""
head = "dtmcs/logs/"
fnames = ["mat_sampling-brp-5994883.out",
          "mat_sampling-herman9-5994914.out",
          "mat_sampling-crowds-5994925.out",
          "mat_sampling-leader-5994880.out",
          "mat_sampling-herman-5994866.out"]

res = [] # length, states, transitinons, precomp, sample
for fname in fnames:
    with open(head+fname) as f:
        content = f.read()
    for match in re.finditer(pat, content):
        line = list(match.groups())
        res.append([e if e else "-1" for e in line])
            #res.append([match.group(i) for i in [2, 3, 4, 5, 6, 7, 8]])

# model, params, length, parse, #states, #transitions,
# precompute, written mats, prob, avg/sample

fname = 'dtmcs/mat_timing.csv'
with open(fname, 'w') as f:
    f.write('name,params,length,parse_time,states,trans,precomp_time,newmat,prob,sample_time\n')
    f.writelines([(','.join([p for p in line]) + '\n') for line in res])

print(f"Written {len(res)} entries to {fname}")