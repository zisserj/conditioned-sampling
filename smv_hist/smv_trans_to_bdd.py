import argparse
import pynusmv


parser = argparse.ArgumentParser("Extracts transitions of SMV system into a BDD.")
parser.add_argument("fname", help="Target SMV file.", type=str)
args = parser.parse_args()
fname = args.fname

#fname="smv_examples/modcounter.m4.smv"

pynusmv.init.init_nusmv()
pynusmv.glob.load_from_file(fname)
pynusmv.glob.compute_model()

fsm = pynusmv.glob.prop_database().master.bddFsm
encoder = fsm.bddEnc
trans_rel = fsm.trans.monolithic

out_file = fname.replace('.smv', '_trans.bdd')

with open(out_file, 'w+') as f:
	encoder.dump(trans_rel, f)
