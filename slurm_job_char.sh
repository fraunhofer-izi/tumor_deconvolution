#! /bin/bash -
#SBATCH --job-name="generate cell type characteristics"
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/genchar.%j.slurmlog
#SBATCH --exclude=ribnode013,ribnode014,ribnode015

source ~dominik.otto/enterEnv.sh # conda environment with python 3.6
python python/genchar.py $*
