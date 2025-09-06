import argparse
import time
import stormpy
import stormpy.simulator

from sparse_mat_sample import ms_str_from, ms_str_any


def sample_trace(sim, length, target):
    tr = []
    obs, reward, labels = sim.restart()
    tr.append(obs)
    for k in range(length):
        if sim.nr_available_actions() > 0:
            obs, reward, labels = sim.step()
            tr.append(obs)
        else:
            return []
    if target in labels:
        return tr
    else:
        return []

def sample_relevant_traces(count, max_r, sim, length, target):
    traces = []
    for i in range(1, max_r+1):
        tr_attempt = sample_trace(sim, length, target)
        if len(tr_attempt) > 0:
            traces.append(tr_attempt)
            if len(traces) == count:
                return i, traces
    return i, traces # type: ignore

if __name__ == "__main__":
    parser = True
    # python sim_sample.py dtmcs/die.pm 8 -repeats 10
    if parser:
        parser = argparse.ArgumentParser("Generates trace samples of system Storm.")
        parser.add_argument("fname", help="Model exported as drn file by storm", type=str)
        parser.add_argument("length", help="Generated trace length", type=int)
        parser.add_argument("-repeats", help="Number of traces to generate ending at target", type=int, default=100)
        parser.add_argument("-tlabel", help="Name of target label matching desired final states",
                            type=str, default='target')
        parser.add_argument('--max-repeats', type=int, help="Maximum number of all traces generated, defaults to repeats*100")
        parser.add_argument('--constants', help="Assignment to undeclared constants", type=str)
        
        args = parser.parse_args()
        filename = args.fname
        path_n = args.length
        repeats = args.repeats
        tlabel = args.tlabel
        max_repeats = args.max_repeats
        const_str = args.constants
    else:
        filename = "dtmcs/nand/nand.pm"
        path_n = 128
        repeats = 100
        tlabel = 'target'
        max_repeats = None
        const_str=None
    if not max_repeats:
        max_repeats = repeats * 100
    
    print(f'Running parameters: fname={filename}, n={path_n}, repeats={repeats}, label={tlabel}, max_repeats={max_repeats}')
    parse_time = time.perf_counter_ns()
    
    prism_program = stormpy.parse_prism_program(filename)
    if const_str:
        desc = stormpy.SymbolicModelDescription(prism_program)
        const_asgn = desc.parse_constant_definitions(const_str)
        prism_program = desc.instantiate_constants(const_asgn).as_prism_program()
    model = stormpy.build_model(prism_program)
    sim = stormpy.simulator.create_simulator(model)
    assert sim
    print(f'Finished creating simulator: {ms_str_from(parse_time)}.')
    
    sim_time = time.perf_counter_ns()
    attempts, res = sample_relevant_traces(repeats, max_repeats, sim, path_n, tlabel)
    final_sim_time = time.perf_counter_ns() - sim_time
    if attempts >= max_repeats:
        print(f"Failed to sample {repeats} conditional traces in {max_repeats} attempts (got {len(res)}) in {ms_str_any(final_sim_time)}.")
    if len(res) > 0:
        print(f'Taken {ms_str_any(final_sim_time/len(res))} per sample')