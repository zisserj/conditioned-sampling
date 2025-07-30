#!/bin/bash

# Usage: ./sample_mat.sh <dir>
# Example: ./sample_mat.sh brp

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_file>"
    echo "Example: $0 brp"
    exit 1
fi

dir="$1"
script="../../sparse_mat_sample.py"

if [ ! -d "$dir" ]; then
    echo "Error: Directory '$dir' not found"
    exit 1
fi

cd "$dir"
for file in *.drn; do
    for length in 8 16 32 64 128; do
        echo "--- $file ---"
        python $script $file $length -repeats 5000
    done
done
# sbatch <<EOT


# #!/bin/bash

# ################################################################################################
# ### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
# ### To ignore, just add another # - like so: ##SBATCH
# ################################################################################################

# #SBATCH --partition main			### specify partition name where to run a job. change only if you have a matching qos!! main: all nodes; gtx1080: 1080 gpu card nodes; rtx2080: 2080 nodes; teslap100: p100 nodes; titanrtx: titan nodes
# #SBATCH --time 0-10:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
# #SBATCH --job-name bench_mats			### name of the job
# #SBATCH --output job-%J.out			### output log for running job - %J for job number
# #SBATCH --gpus=0				### number of GPUs, allocating more than 1 requires IT team's permission. Example to request 3090 gpu: #SBATCH --gpus=rtx_3090:1

# # Note: the following 4 lines are commented out
# ##SBATCH --mail-user=user@post.bgu.ac.il	### user's email for sending job status messages
# ##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
# ##SBATCH --mem=24G				### ammount of RAM memory, allocating more than 60G requires IT team's permission

# #SBATCH -o "outFile"$1".txt"
# #SBATCH -e "errFile"$1".txt"

# ################  Following lines will be executed by a compute node    #######################

# ### Print some data to output file ###
# echo `date`
# echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
# echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

# ### Start your code below ####
# module load anaconda				### load anaconda module (must be present when working with conda environments)
# source activate conditional				### activate a conda environment, replace my_env with your conda environment

# cd ~/conditioned/graph-sampling/dtmcs
# for dir in brp crowds egl herman leader_sync nand; do
#     for length in 8 16 32 64 128; do
#         pushd $dir
#         bash ../bench_mat.sh $dir
#         echo done $file
#         popd
# done

# ### run my_code.py python file and send my_arg1 argument to it
