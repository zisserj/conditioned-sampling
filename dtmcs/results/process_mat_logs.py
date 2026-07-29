import re

entry_pat = r"--- \w+\/([a-zA-Z]+)_?([\S_]+).drn - (\d+) ---\n"

output_pat = r'''[\w ]+input: ([\d.]+)ms.\
[\w ]+states: ([\d.]+)\
[\w ]+transitions: ([\d.]+)\
(?:\s?Found [\w\W]+?\n?)?(?:[\w ]+functions: ([\d.]+)ms.)?
(?:[\w ]+ mats: G\(\d+\), T\((\d+)\)\n)?[\w ]+ probability is ([\d.+-e]+)\
(?:[#\d\-\s]*Taken ([\d.]+)ms per sample|No matching traces)'''



content = ""
dir = "dtmcs/features/"
fnames = '''mat_sampling-herman_mat-19611607.out
mat_sampling-herman_mat-19621651.out'''.split()


res = [] # length, states, transitinons, precomp, sample
for fname in fnames:
    with open(dir+fname) as f:
        content = f.read()
    seq = re.split(entry_pat, content) # [preamble, model, params, length, content, model,...]
    for i in range(1, len(seq), 4):
        name, params, length, output_content = seq[i:i+4]
        output_type = "ok"
        stats = ["-1" for _ in range(7)]
        if "OverflowError" in output_content:
            output_type = "overflow"
        elif "Killed" in output_content:
            output_type = "Killed"
        elif "OOM Killed" in output_content:
            output_type = "oom"
        elif "ValueError: a cannot be empty unless no samples are taken" in output_content:
            output_type = "numeric"
        else:        
            match = re.search(output_pat, output_content, re.M)
            if match:
                print(match.groups())
                # parsetime, #states, #transitions,
                # precompute, written mats, prob, avg/sample
                stats = [e if e else "-1" for e in match.groups()]
            else:
                print(f"Issue processing {name}-{params} ({length}): {output_content}")
        entry = [name, params, length] + [output_type] + stats
        res.append(','.join(entry))
# model, params, length, parsetime, #states, #transitions,
# precompute, written mats, prob, avg/sample

fname = dir+'mat_timing.csv'
with open(fname, 'w') as f:
    f.write('name,params,length,output_type,parse_time,states,trans,precomp_time,newmat,prob,sample_time\n')
    f.write('\n'.join(res))

print(f"Written {len(res)} entries to {fname}")