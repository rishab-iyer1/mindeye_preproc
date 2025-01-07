#!/bin/bash

# Name of job?
#SBATCH --job-name=GLMsingle-multi-23-ses-03

# Set partition
#SBATCH --partition=all

# How long is job?
#SBATCH -t 05:00:00

# How much memory to allocate (in MB)?
#SBATCH --cpus-per-task=14 --mem-per-cpu=24000

#SBATCH -e slurms/%j.err        # first create a "slurms" folder in current directory to store logs
#SBATCH -o slurms/%j.out

#SBATCH --mail-user=rsiyer@princeton.edu
#SBATCH --mail-type=BEGIN,END

module purge
source ~/fmri/bin/activate
cd /usr/people/ri4541/rtmindeye/code/analysis

jupyter nbconvert --to script GLMsingle-multisession-2-3.ipynb && ipython GLMsingle-multisession-2-3.py
