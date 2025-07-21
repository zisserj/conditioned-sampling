import pynusmv
pynusmv.init.init_nusmv()
# pynusmv.glob.load_from_file("/home/jules/benchmark/examples/example_cmu/counter.smv")
# pynusmv.glob.load_from_file("/home/jules/benchmark/examples/bmc_tutorial/bmc_tutorial.smv")
pynusmv.glob.load_from_file("smv_examples/modcounter.m4.smv")

pynusmv.glob.compute_model()

fsm = pynusmv.glob.prop_database().master.bddFsm
encoder = fsm.bddEnc

trans_rel = fsm.trans.monolithic
ins =  pynusmv.dd.Inputs.from_bdd(trans_rel, fsm)
# print(ins.get_str_values())



with open('mod_trans.bdd', 'w+') as f:
	encoder.dump(trans_rel, f)

with open('mod_reachable.bdd', 'w+') as f:
	encoder.dump(fsm.init, f)

newfsm = pynusmv.fsm.BddFsm(trans_rel._ptr)
