#!/bin/bash
#SBATCH --account p31225  ## <-- EDIT THIS TO BE YOUR ALLOCATION 
#SBATCH --partition short  ## <-- EDIT THIS TO BE YOUR QUEUE NAME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=sample_job
#SBATCH --output=outlog
#SBATCH --error=errlog

module purge all                ## Unload existing modules
module load python  		## Load necessary modules (software, libraries)

bash /projects/intro/whereami.sh         ## Run the program
python /projects/intro/helloworld.py     ## Run the program
