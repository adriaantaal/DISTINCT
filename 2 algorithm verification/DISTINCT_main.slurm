#!/bin/bash
#SBATCH --job-name=DSonly    		# Job name
#SBATCH --mail-type=FAIL        		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=at3111@columbia.edu     # Where to send mail
#SBATCH --ntasks=1                  		# Run on a single CPU (matlab can only use ntasks=1)
#SBATCH --mem=200gb              			# Job memory request
#SBATCH --time=168:00:00         			# Time limit hrs:min:sec
#SBATCH --output="matlab_distinct.log"   	# Standard output and error log
#SBATCH --qos=ks_distributedcomputing 
#SBATCH --partition=charmander

source ~/.bashrc;

pwd; hostname; date
echo "Running plot script"
matlab -nosplash -noFigureWindows -nodesktop -r "try; run('/proj2/adriaantaal/96 Paper DISTINCT/2 algorithm verification/DISTINCT_main_DSonly'); catch; end; quit;"