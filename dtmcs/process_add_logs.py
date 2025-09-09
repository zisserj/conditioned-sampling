import re

logs_param = r"--- \w+\/([a-zA-Z]+)_?([\S_]+)\.drdd - (\d+) ---\
([\s\S]*?)(?=--- \w+\/[a-zA-Z]+_?[\S_]+\.drdd - \d+ ---|$)"


output_pat = r"[\w ]+input: ([\d.]+)ms.\
[\w ]+state: ([\d.]+)\
[\w ]+ADD: ([\d.]+) nodes\
[\w ]+functions: ([\d.]+)ms.\
(?:[#\d\-\s]*Taken ([\d.]+)ms per sample|No matching traces)"


content = ""
head = "dtmcs/logs/"
fnames = ["add_sampling-crowds-6482564.out",
          "add_sampling-brp-6482563.out",
          "add_sampling-egl-6482658.out",
          "add_sampling-herman-6482661.out",
          "add_sampling-leader_sync-6482666.out",
          "add_sampling-nand-6482668.out"]
for fname in fnames:
    with open(head+fname) as f:
        content += f.read()

res = [] # model, params, length, output
for fname in fnames:
    with open(head+fname) as f:
        content = f.read()
    for match in re.finditer(logs_param, content):
        name, params, length, output_content = match.groups()
        output_type = "ok"
        parsetime = varnum = addsize = precomptime = tracetime = '0'
        if "Segmentation fault" in output_content:
            output_type = "segfault"
        elif "CANCELLED" in output_content:
            output_type = "timeout"
        elif "CUDD appears to have run out of memory." in output_content:
            output_type = "mem"
        else:
            try:
                stats = re.findall(output_pat, output_content)[0]
                parsetime, varnum, addsize, precomptime, tracetime = stats
                if not tracetime:
                    tracetime = "-1"
            except:
                print("Issue processing", name, params, length, " : ", output_content)
        res.append(','.join([name, params, length, output_type,
                            parsetime, varnum, addsize, precomptime, tracetime]))

fname = 'dtmcs/add_timing.csv'
with open(fname, 'w') as f:
    f.write("name,params,length,output_type,parsetime,vars_per_state,add_size,precomp_time,trace_time\n")
    f.write('\n'.join(res))

print(f"Written {len(res)} entries to {fname}")