import re

pat = r"--- ([a-zA-Z]+)_?([\S_]*).pm (?:--constants )?([\w=,]*) ?- (\d+) ---\
[\s\S]*?repeats=(\d+)[\s\S]*?\
[\w ]+simulator: ([\d.]+)ms\.\
(?:[#\d\-\s]*Taken ([\d.]+)ms per sample|[\w\d\s]+ \(got (\d+)\) in ([\d.]+)ms.)"


content = ""
head = "dtmcs/logs/"
fnames = ["sim_sampling-crowds-6483701.out",
          "sim_sampling-leader-6483709.out",
          "sim_sampling-egl-6483703.out",
          "sim_sampling-herman-6483711.out",
          "sim_sampling-nand-6483706.out",
          "sim_sampling-brp-6483695.out"]
for fname in fnames:
    with open(head+fname) as f:
        content += f.read()

with open('test.txt', 'w') as f:
    f.write(content)
    
res = [] # model, params, length, output
for fname in fnames:
    with open(head+fname) as f:
        content = f.read()
    for match in re.finditer(pat, content):
        print(match.groups())
        name, params_a, params_b, length, repeats, setup_time, sample_time, suc_count, total_time  = match.groups()
        params = params_a if params_a else params_b.replace(',', ' ')
        if sample_time:
            suc_count = repeats
            total_time = "-1"
        else:
            sample_time = "-1"
        res.append(','.join([name, params, length, setup_time,
                            sample_time, suc_count, total_time]))

print(res)
fname = 'dtmcs/sim_timing.csv'
with open(fname, 'w') as f:
    f.write("name,params,length,setup_time,sample_time,suc_count,total_time\n")
    f.write('\n'.join(res))

print(f"Written {len(res)} entries to {fname}")