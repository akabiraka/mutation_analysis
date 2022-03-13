#!/usr/bin/sh

## this must be run from directory where run.py exists.
## --workdir is not used in this file.

#SBATCH --job-name=mutant
#SBATCH --output=/scratch/akabir4/mutation_analysis/outputs/argo_logs/pssm-%N-%j.output
#SBATCH --error=/scratch/akabir4/mutation_analysis/outputs/argo_logs/pssm-%N-%j.error
#SBATCH --mail-user=<akabir4@gmu.edu>
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --partition=all-LoPri
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000MB

#SBATCH --array=1315-1318
python generators/distributed_pssm.py
