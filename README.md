## *MCCS: Markov Chain Conditioned Sampling*

This repository features three realisations of an algorithm which generates sample traces of specified length of an input Markov chain, and meet the provided initial and end conditions.
It was implemented as part of my MSc thesis research.

The realisations are based on different data structures: Algebraic Decision Diagrams (ADDs), Binary Decision Diagrams (BDDs), and sparse tensors. All the files related to implementation are in this main directory, and those used for the evaluation are in `dtmcs`.

### Pipeline usage 

<img src="cond_sampling_pipeline.png" width="75%"/>

Default usage meant for [Prism](https://www.prismmodelchecker.org) DTMC models, like those found in the `examples` folder. The set of initial states is assumed from the model, targeted states must be assigned a label.

**Model processing** requires [Storm](https://www.stormchecker.org/index.html), and is performed by the commands
```Bash
storm --prism dice.pm --buildfull --prismcompat --engine sparse --exportbuild dice.drn
storm --prism dice.pm --buildfull --prismcompat --engine dd --exportbuild dice.drdd
```
Note the usage of different `--engine` options depending on desired output type: `sparse` for `drn` (matrices), and `dd` for `drdd` (ADD/BDD). Constants assignments can be passed as parameters, for example for `N` and `MAX` use `--constants N=16,MAX=2`.

Some processed files have been provided in the `example` folder for convenience.

**Script usage** of the sampling modules `x_sample.py` is largely the same across implementations, and argument help (`-h`) is provided. These scripts call the parsing functions within `x_to_y.py` automatically, there is no need to use them explicitly.
`sim_sample.py`, unlike the others, performs naive rejection sampling as opposed to our conditioned algorithm and takes a *prism* file as input Note it requires stormpy.

```
python sparse_mat_sample.py --help
usage: Generates conditional samples of system via sparse matrices. [-h] [-repeats REPEATS] [-tlabel TLABEL]
                                                                    [-output OUTPUT] [--alg4] [--store]
                                                                    fname length
positional arguments:
  fname             Model exported as drn file by storm
  length            Generated trace length

options:
  -h, --help        show this help message and exit
  -repeats REPEATS  Number of traces to generate
  -tlabel TLABEL    Name of target label matching desired final states
  -output OUTPUT    File destination for generated traces
  --alg4            Use general alg for powers of two instead of baseline
  --store           Store / try loading existing mats
```

And example usage:
```
python add_sample.py robot.drdd 32 -repeats 100 -output robot.drdd.out
```


### Direct usage

The algorithms can be invoked directly by importing the relevant modules. This is easiest to do with the sparse tensor implementation, as they are easiest to manipulate directly. While the exact functions used depend on the implementation, the general workflow should be as follows:
1. Define transition, initial and target set functions.
2. Precompute the multi-step transition functions via the `compute_power` functions.
3. Generate samples using the information from the previous step via `draw_sample` or `generate_many_traces`.

The `__main__` functions of each module can be used as a reference, and otherwise feel free to reach out for questions. Matching common Python conventions, the underscored functions are considered internal and should not be invoked directly.

