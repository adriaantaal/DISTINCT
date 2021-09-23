#!/bin/bash
#SBATCH --job-name=DSONLY		    		# Job name
#SBATCH --mail-type=FAIL        		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=at3111@columbia.edu     # Where to send mail
#SBATCH --ntasks=1                  		# Run on a single CPU (matlab can only use ntasks=1)
#SBATCH --cpus-per-task=1             		# Run on one CPUs
#SBATCH --mem=200gb              			# Job memory request. max 200gb per user
#SBATCH --time=168:00:00         			# Time limit hrs:min:sec
#SBATCH --output="matlab_distinct.log"   	# Standard output and error log
#SBATCH --qos=ks_lowprio 
#SBATCH --partition=pikachu

NUMBATCHFILE=${1} # Note that this line gets replaced, should always be line 13 in the script

SLURM_OUTPUT_FILE=results/slurm_user_${NUMBATCHFILE}.log
SLURM_WORKING_DIR=$(pwd)
SLURM_HOST=$(hostname)
SLURM_DATE=$(date)

source ~/.bashrc
echo -e "Working directory is $SLURM_WORKING_DIR\n" > $SLURM_OUTPUT_FILE
echo -e "Hostname is $SLURM_HOST\n" >> $SLURM_OUTPUT_FILE
echo -e "Date is $SLURM_DATE\n" >> $SLURM_OUTPUT_FILE
echo -e "Running on one CPU \n" >> $SLURM_OUTPUT_FILE
echo -e "Batch file $NUMBATCHFILE\n" >> $SLURM_OUTPUT_FILE

matlab -nosplash -noFigureWindows -nodesktop -r "try; numBatchFile=${NUMBATCHFILE}; run('/proj2/adriaantaal/96 Paper DISTINCT/2 algorithm verification/DISTINCT_main_DSonly'); catch; end; quit;"

exit
