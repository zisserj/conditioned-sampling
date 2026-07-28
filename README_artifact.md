

The artifact submitted is the repository containing the source code, scripts and proccessed logs used for the evaluation portion in our paper.
We have included a Docker file which constructs the work environment used to run said scripts (for transparency sake: the Docker image does download a live version of the Github repository and not the frozen/submitted artifact, but the differences, if any, are likely minor).

Experiments were ran using 32G RAM on a single core, but the the minimal examples should run fine if the PC can manage Docker.

The container can be run (after installing Docker) by using the following commands:
```
docker pull movesrwth/storm:stable
cd <directory of the artifact containing Dockerfile>
docker image build -t mccs:latest
docker run -it mccs:latest /bin/bash
```

The directory `dtmcs/results/fmcad` contains the aggregated results featured in the paper.
`README.md` contains further description of the repository content and usage. Note the BDD implementation is irrelevant to this paper.


An example of commands which generate samples of length 12 for a single model with both backends: uses Storm to export the models as data structures, then the Python scripts for sampling.
```
storm --prism dtmcs/brp/brp.pm --buildfull --prismcompat --engine sparse --constants N=16,MAX=2 --exportbuild dtmcs/brp/brp_N_16_MAX_2.drn
storm --prism dtmcs/brp/brp.pm --buildfull --prismcompat --engine dd --constants N=16,MAX=2 --exportbuild dtmcs/brp/brp_N_16_MAX_2.drdd
python sparse_mat_sample.py dtmcs/brp/brp_N_16_MAX_2.drn 12 -output dtmcs/brp/sparce_traces_12.txt
python add_sample.py dtmcs/brp/brp_N_16_MAX_2.drdd 12 -output dtmcs/brp/add_traces_12.txt
```
And the output:
```
Storm 1.12.0

Date: Wed May 13 09:50:46 2026
Command line arguments: --prism dtmcs/brp/brp.pm --buildfull --prismcompat --engine sparse --constants 'N=16,MAX=2' --exportbuild dtmcs/brp/brp_N_16_MAX_2.drn
Current working directory: /opt/MCCS

Time for model input parsing: 0.007s.

Time for model construction: 0.034s.

--------------------------------------------------------------
Model type:     DTMC (sparse)
States:         677
Transitions:    867
Reward Models:  none
State Labels:   3 labels
   * init -> 1 item(s)
   * deadlock -> 35 item(s)
   * target -> 32 item(s)
Choice Labels:  none
--------------------------------------------------------------

Exporting model to 'dtmcs/brp/brp_N_16_MAX_2.drn'.
Write to file dtmcs/brp/brp_N_16_MAX_2.drn.
Time for model export: 0.001s.

Storm 1.12.0

Date: Wed May 13 09:50:46 2026
Command line arguments: --prism dtmcs/brp/brp.pm --buildfull --prismcompat --engine dd --constants 'N=16,MAX=2' --exportbuild dtmcs/brp/brp_N_16_MAX_2.drdd
Current working directory: /opt/MCCS

Time for model input parsing: 0.002s.

Using Sylvan with 16 parallel threads.
Time for model construction: 1.207s.

--------------------------------------------------------------
Model type:     DTMC (symbolic)
States:         677 (576 nodes)
Transitions:    867 (2182 nodes)
Reward Models:  none
Variables:      rows: 18 meta variables (32 DD variables), columns: 18 meta variables (32 DD variables)
Labels:         3
   * deadlock -> 35 state(s) (106 nodes)
   * init -> 1 state(s) (33 nodes)
   * target
--------------------------------------------------------------
Write to file dtmcs/brp/brp_N_16_MAX_2.drdd.
Running parameters: fname=dtmcs/brp/brp_N_16_MAX_2.drn, n=12, repeats=1000, label=target, store=False, output=dtmcs/brp/sparce_traces_12.txt
Finished parsing input: 8.861257ms.
Number of states: 677
Number of transitions: 867
Finished precomputing functions: 114.211699ms.
Property probability is 1.176e-05
Taken 1.763919ms per sample
1000 traces written to dtmcs/brp/sparce_traces_12.txt.
Running parameters: fname=dtmcs/brp/brp_N_16_MAX_2.drdd, n=12, repeats=1000, label=target, store=False, output=dtmcs/brp/add_traces_12.txt
Finished parsing input: 186.836544ms.
Number of variables per state: 32
Size of ADD: 1391 nodes
Finished precomputing functions: 780.977641ms.
Taken 4.103047ms per sample
1000 traces written to dtmcs/brp/add_traces_12.txt.
```

