#!/bin/bash

################################################################################################
### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --time 0-20:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name herman-a ### name of the job
#SBATCH --output logs/add_sampling-%x-%J.out			### output log for running job - %J for job number
#SBATCH --mem=32G				### ammount of RAM memory, allocating more than 60G requires IT team's permission
#SBATCH --constraint=cpu128

################  Following lines will be executed by a compute node    #######################

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate conditioned_sampling				### activate a conda environment, replace my_env with your conda environment

script="$HOME/conditioned/graph-sampling/add_sample.py"
dir=herman;

for dir in $dir; do # brp crowds  egl herman leader_sync nand; do
  for file in "$dir"/*.drdd; do
    fdir=${file%.*}
    for length in 8 12 16 24 32 48 64 96 128; do
        echo "--- $file - $length ---"
        python $script "$file" $length -repeats 500 -tlabel stable
    done
  done
done
