#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -p shared,stats
#SBATCH --mem=90000
#SBATCH -J job_simple
#SBATCH -o job_simple_%j_%N.out
#SBATCH -e job_simple_%j_%N.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=debug.pie@gmail.com

module load Anaconda3/5.0.1-fasrc02
source activate my_root

module load GCC/8.2.0-2.31.1
module load R/3.6.3-fasrc01
# module load Python/3.7.2

Rscript --vanilla  ../others_timer.R --method ${1} --cnm ${2} --ncores 8
