#!/bin/bash

# Specify the number of processes to create
START_FILE_NUM=$1
END_FILE_NUM=$2

for ((i=${START_FILE_NUM}; i <= ${END_FILE_NUM} ; i++));
do
	echo "Launching sbatch process ${i}"
	awk -v my_var="${i}" 'NR==13 {$0="NUMBATCHFILE=" my_var} 1' DISTINCT_singleinstance.slurm > gen_slurm_script_${1}.sh
	chmod 744 gen_slurm_script_${1}.sh
	sbatch gen_slurm_script_${1}.sh
	rm gen_slurm_script_${1}.sh
done
