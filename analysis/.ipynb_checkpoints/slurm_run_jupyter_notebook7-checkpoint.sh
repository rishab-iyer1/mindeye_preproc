#!/bin/bash

# Name of job?
#SBATCH --job-name=GLMsingle

# Set array to be your subject number
#SBATCH --array=001


# Where to output log files?
#SBATCH --output='./logs/assessGLMsingle_sub%a_%A.log'

# Set partition
#SBATCH --partition=all

# How long is job?
#SBATCH -t 12:00:00

# How much memory to allocate (in MB)?
#SBATCH --cpus-per-task=16 --mem-per-cpu=14000

#SBATCH --mail-user=eam7@princeton.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Print job submission info
echo "Slurm job ID: " $SLURM_JOB_ID
date

# Set subject ID based on array index
printf -v subject "%03d" $SLURM_ARRAY_TASK_ID

# Activate conda environment
module purge #make sure no modules are loaded
source /usr/people/eam7/.bashrc
conda activate mindeye

# RUN NOTEBOOK(S)

# notebook to run:
ntbk=07
notebook=./GLMsingle_mindeye_sub-${subject}.ipynb

# output directory:
outDir=./notebook-glmsingle-output

echo "Running notebook $ntbk for sub-$subject"
# first notebook:
jupyter nbconvert \
  --ExecutePreprocessor.allow_errors=True \
  --ExecutePreprocessor.timeout=-1 \
  --FilesWriter.build_directory=$outDir \
  --to html \
  --execute $notebook

# after running notebook, move .html from temp output directory to permanent location
temp=./notebook-glmsingle-output/*${subject}.html
outDir=../../data/bids/derivatives/glmsingle/GLMsingle_mindeye_sub-${subject}.html

# mv from temp directory to permanent directory:
mv $temp $outDir
echo "Finished running notebook $ntbk for sub-$subject"


# # NEXT NOTEBOOK
# ntbk=05 #notebook number for output directory

# echo "Running notebook $ntbk for sub-$subject"

# ./run_jupyter_notebook.sh $subject

# # after running notebook, move .html from temp output directory to permanent location
# # mv from temp output directory:
# temp=./notebook-${ntbk}-output/*${subject}.html

# # mv to permanent directory:
# outDir=../../data/mainanalysis/notebook-${ntbk}-output/${ntbk}-PSanalysis-v4masks-thresh75_sub-${subject}.html

# mv $temp $outDir

# echo "Finished running notebook $ntbk for sub-$subject"

date
