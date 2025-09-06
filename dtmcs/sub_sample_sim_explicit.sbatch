#!/bin/bash

# Usage: sbatch sub_sample_sim.sbatch <input_file>
# Example: sbtach sub_sample_sim.sbatch brp/models


################################################################################################
### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main                        ### specify partition name where to run a job. change only if you have a matching qos!! main: all nodes; gtx1080: 1080 gpu card nodes; rtx2080: 2080 nodes; teslap100: p100 nodes; titanrtx: titan nodes
#SBATCH --time 0-20:30:00                       ### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name herman ### name of the job
#SBATCH --output logs/sim_sampling-%x-%J.out                    ### output log for running job - %J for job number
#SBATCH --gpus=gtx_1080:1
# Note: the following 4 lines are commented out
##SBATCH --mail-user=user@post.bgu.ac.il        ### user's email for sending job status messages
##SBATCH --mail-type=ALL                        ### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=32G                               ### ammount of RAM memory, allocating more than 60G requires IT team's permission


################  Following lines will be executed by a compute node    #######################

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda                            ### load anaconda module (must be present when working with conda environments)
source activate conditioned_sampling                            ### activate a conda environment, replace my_env with your conda environment


if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_file>"
    echo "Example: $0 brp/models"
    exit 1
fi

input_file="$1"
script="$HOME/conditioned/graph-sampling/sim_sample.py"
label=stable;

if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found"
    exit 1
fi

while IFS="" read -r line || [ -n "$line" ]
do
    # Skip empty lines and comment lines (starting with #)
    if [ -n "$line" ] && [ "${line#\#}" = "$line" ]; then
        # Extract first and third words from the line
        set -- $line
	dir_name=$(echo $1 | sed 's/[0-9].*//')
	dir_name="${dir_name//.pm/}"
	for length in 8 16 32 64 128; do
            echo "--- $line - $length ---"
            python $script "$dir_name"/$line $length -repeats 500 -tlabel "$label"
        done

        fi
done < "$input_file"
