#! /bin/bash -
#SBATCH --job-name="feature selection"
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/featureselection.%j.slurmlog

source /etc/profile.d/modules.sh
source ~dominik.otto/enterEnv.sh # conda environment with python 3.6
module load R/3.5.2-0

python python/unify_expression_HUGO_greed.py $*
