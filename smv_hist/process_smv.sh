#!/bin/bash
set -e

smv_fname="$1"
bdd_fname="${smv_fname/.smv/_trans.bdd}"
mat_fname="${bdd_fname/.bdd/.smv}"

nusmv_py="/home/jules/.env/nusmv/bin/python3.6"
reg_py=".venv/bin/python"

echo ${nusmv_py} smv_trans_to_bdd.py ${smv_fname}

${nusmv_py} smv_trans_to_bdd.py ${smv_fname}
${reg_py} bdd_to_mat.py ${bdd_fname}
${reg_py} mat_sample.py ${mat_fname}