#!/bin/bash
#SBATCH -A p31225						# allocation
#SBATCH -p normal						# queue partition assignment, 'normal' or 'short'
#SBATCH -t 08:00:00						# duration of job
#SBATCH -N 1							# number of computers needed
#SBATCH --array=1-180 						# IMPORTANT: number of parallel jobs, per the lines in the csv file
#SBATCH --ntasks-per-node=1 					# how many cpus do you need on each computer
#SBATCH --mem-per-cpu=900M 					# how much RAM needed per CPU (affects FairShare score so be careful)
#SBATCH --job-name="SIMRUN_\${SLURM_ARRAY_TASK_ID}"	 	# use the task id in the name of the job
#SBATCH --output=SIMRUN.%A_%a.out 				# jobid (A) and job index (a) to name the log file
#SBATCH --mail-type=BEGIN,END					# for email alertss, one of BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --mail-user=eugenewickett@u.northwestern.edu  		# your email

module purge all
# module load python/anaconda
module load python-anaconda3/2019.10
# module load python/ActivePython-3.2
#source activate slurm-py37-test

# home/eow3176/... (?)
cd SFP_Sim/Falsification_Simulation_Model
python Falsification_Simulator_QUEST.py $(sed -n ${SLURM_ARRAY_TASK_ID}p parameterFile.csv | sed 's/,/ /g')

# make this file executable and then run from the command line
# chmod u+x [script name]
# ./[script name]
