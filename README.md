# Princeton 3T MindEye Data Processing Pipeline
## Pre-processing

1. data transfer
    - ```ssh ri4541@scotty.princeton.edu```
    - Create a new folder called "data" and subfolders called “dicom” and “work” 
        - this way you don't have to modify globals.sh; it looks for a folder named "data"
        - later, rename with BIDS name e.g. data_sub-003_ses-01
    - transfer dicoms from jukebox to the newly created data folder
        - example usage: ```rsync -aP /jukebox/dicom/conquest/Prisma-MSTZ400D/NormaL/2024/003_ses01_rtmindeye-1213-1401/ /jukebox/norman/rsiyer/rtmindeye/data/dicom```
    - Rename the dcm folder to the dicom's name, e.g. 03_ses01_rtmindeye-1213-1401


2. setup
    - ```cd /Volumes/norman/rsiyer/rtmindeye/code```
    - ```module load fsl``` needed for PyDeface
    - ```source ~/fmri/bin/activate``` PyDeface should be pip installed 

1. ```preprocessing/step1_preproc.sh```
    - ```cd preprocessing```
    - example usage: ```./step1_preproc.sh 003 01 003_ses01_rtmindeye-1213-1401```
        - uses heudiconv to convert DICOMs to BIDS format
        - uses PyDeface to deface images
      
2. ```preprocessing/step2_preproc.sh```
    - example usage: ```./step2_preproc.sh 003 01```
    - delete scouts and duplicate runs from the .tsv dile to match the files that are deleted by the step2 script
    - modify IntendedFor field in the fieldmap JSON files to be valid JSON format
    - current issues (TODO)
        - IntendedFor field in fmap folder gets added repeatedly if it already exists causing JSON errors
        - .tsv file doesn't update with deleted scout files
      
3. [BIDS Validator](https://bids-standard.github.io/bids-validator/)
    - verify manually that the data is BIDS-compatible
        - address any red (errors); can ignore the yellow (warnings) because fMRIPrep will still work
      
4. ```preprocessing/run_mriqc.sh```
    - example usage: ```./run_mriqc.sh 003```
    - quality control checks, verify outputs manually to look for outliers 
        - absolute values don't matter as much
      
5. ```preprocessing/run_mriqc_group.sh```
    - example usage: ```./run_mriqc_group.sh```
    - quality control checks, verify outputs manually to look for outliers
        - absolute values don't matter as much

6. ```preprocessing/run_fmriprep.sh```
    - runs fMRIPrep for all sessions for an individual participant


## GLMsingle
1. setup
    - create data/bids/derivatives/masks/sub-003 if it doesn't exist
    - create data/design/csv if it doesn't exist
    - populate csv folder with the run-by-run csv files from GitHub
        - for example: https://github.com/PrincetonCompMemLab/real_time_mindEye2/tree/paul/psychopy_task/conditions_files/participant3_run0_sess1.csv
    - also populate csv folder with design file from GitHub
        - for example: https://github.com/PrincetonCompMemLab/real_time_mindEye2/blob/paul/psychopy_task/data/3_1_rtmindeye_2024-12-13_14h18.09.698.csv
    - rename the design file from 3_1_rtmindeye_2024-12-13_14h18.09.698.csv to sub-003_ses-01.csv

2. ```analysis/nsdgeneral_to_epi.sh```
    - create an NSDgeneral mask for a subject containing primarily visual cortex voxels
      
3. ```analysis/GLMsingle.ipynb```
    - get model accuracy (r^2)
    - get single-trial betas
    - get best-fitting HRF
    - get fractional ridge regression regularization level
    - verify reliability and quality control plots manually

## MindEye
TODO: describe where to copy the GLMsingle outputs and masks to the workstation and how to run mindeye