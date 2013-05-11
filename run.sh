#!/bin/bash
#SBATCH --job-name=a1
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -n 1

date
srun ./feaprog
date
