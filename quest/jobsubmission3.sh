#!/bin/bash

# this next line reads tab-delimited (\t) and expects to get two values per line, 
# saved as P1 and P2
while IFS=$'\t' read P1 P2 P3 P4
do
    JOB=`sbatch << EOJ
#!/bin/bash
#SBATCH --account=p31225  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal  ### PARTITION (buyin, short, normal, w10001, etc)
#SBATCH --array=0 ## number of jobs to run "in parallel"
#SBATCH --nodes=1  ## how many computers do you need
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on each computer
#SBATCH --time=03:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=900M ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than y$
#SBATCH --job-name="sample_job_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=sample_job.%A_%a.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, et$
#SBATCH --mail-user=eugenewickett@u.northwestern.edu  ## your email
module purge all
# module load python/anaconda
module load python-anaconda3/2019.10
# module load python/ActivePython-3.2
#source activate slurm-py37-test

# Here goes the application-specific line.  Example here is using 
# named command line arguments for a Python script, but use ${P1} and ${P2} a
# as appropriate for your application
python python_finder/Falsification_Simulator.py --diag_sens=${P1} --diag_spec=${P2} --mult_days=${P3} --sim_days=${P4}
EOJ
`

done < python_finder/parameterFile.txt
exit   

# make this file executable and then run from the command line
# chmod u+x submit.sh
# ./submit.sh
