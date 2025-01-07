#!/bin/bash

# Name of job?
#SBATCH --job-name=sub-003

# Set partition
#SBATCH --partition=all

# How long is job?
#SBATCH -t 20:00:00

# How much memory to allocate (in MB)?
#SBATCH --cpus-per-task=10 --mem-per-cpu=16000

#SBATCH -e slurms/%j.err        # first create a "slurms" folder in current directory to store logs
#SBATCH -o slurms/%j.out

#SBATCH --mail-user=rsiyer@princeton.edu
#SBATCH --mail-type=END

module purge
source ~/mindeye/bin/activate
cd /usr/people/ri4541/rtmindeye/code/analysis

jupyter nbconvert --to script GLMsingle.ipynb && ipython GLMsingle.py
