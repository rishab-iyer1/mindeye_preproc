# Princeton 3T MindEye Data Processing Pipeline
## Pre-processing

1. data transfer
    - ```ssh ri4541@scotty.princeton.edu```

    - Create a new folder called "data" and subfolders called “dicom” and “work” 
        - this way you don't have to modify globals.sh; it looks for a folder named "data"
        - later, rename with BIDS name e.g. data_sub-003_ses-01
    - ```tmux new -s preproc``` use [screen](https://gist.github.com/jctosta/af918e1618682638aa82) or [tmux](https://tmuxcheatsheet.com/) to keep commands running in the background
        - do all the preprocessing in this environment so it's easier to trace back; it stays up until you delete the session
    - transfer dicoms from jukebox to the newly created data folder
        - create a directory inside dicom called e.g. 003_ses01_rtmindeye-1213-1401
        - example usage: ```rsync -aP /jukebox/dicom/conquest/Prisma-MSTZ400D/NormaL/2024/003_ses01_rtmindeye-1213-1401 /jukebox/norman/rsiyer/rtmindeye/data/dicom/```


2. setup
    - ```cd /jukebox/norman/rsiyer/rtmindeye/code/```
    - ```module load fsl``` needed for PyDeface
    - ```source ~/fmri/bin/activate``` 
    - PyDeface should be pip installed prior to continuing

3. ```preprocessing/step1_preproc.sh```
    - ```cd preprocessing```
    - example usage: ```./step1_preproc.sh 003 01 003_ses01_rtmindeye-1213-1401```
        - uses heudiconv to convert DICOMs to BIDS format
        - uses PyDeface to deface images
      
4. ```preprocessing/step2_preproc.sh```
    - example usage: ```./step2_preproc.sh 003 01```
    - delete scouts and duplicate runs from the .tsv dile to match the files that are deleted by the step2 script
    - modify IntendedFor field in the fieldmap JSON files to be valid JSON format
    - current issues (TODO)
        - IntendedFor field in fmap folder gets added repeatedly if it already exists causing JSON errors
        - .tsv file doesn't update with deleted scout files
     
5. [BIDS Validator](https://bids-standard.github.io/bids-validator/)
    - verify manually that the data is BIDS-compatible
        - address any red (errors); can ignore the yellow (warnings) because fMRIPrep will still work

Can easily run the next three commands in sequence with ```./run_fmriprep.sh 003 && ./run_mriqc.sh 003 && ./run_mriqc_group.sh```      


6. ```preprocessing/run_mriqc.sh```
    - example usage: ```./run_mriqc.sh 003```
    - quality control checks, verify outputs manually to look for outliers 
        - absolute values don't matter as much
      
7. ```preprocessing/run_mriqc_group.sh```
    - example usage: ```./run_mriqc_group.sh```
    - quality control checks, verify outputs manually to look for outliers
        - absolute values don't matter as much

8. ```preprocessing/run_fmriprep.sh```
    - example usage: ```./run_fmriprep.sh 003```
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

2. Rename the data folder e.g. data_sub-003_ses-01

3. ```analysis/nsdgeneral_to_epi.sh```
    - ```cd ../analysis/```
    - example usage: ```./nsdgeneral_to_epi.sh 005 ses-04 D _task-D --data-folder=data_sub-005_ses-04```
    - create an NSDgeneral mask for a subject containing primarily visual cortex voxels
      
4. ```analysis/GLMsingle.ipynb```
    - if submitting a SLURM job, edit ```run_glmsingle.sh``` and ```GLMsingle.ipynb``` with the desired paths and settings
    - if running directly (say on Scotty's login node):
        - ```source ~/mindeye/bin/activate```
        - use ```run_glmsingle.py``` to easily set various environment variables as command line arguments
        - ```run_glmsingle.py --help``` for a list of all accepted arguments
        - example usage: ```python ~/rtmindeye/code/analysis/run_glmsingle.py data_sub-005_ses-04 glmsingle_ses-04_task-D sub-005 ses-04 D```
    - get model accuracy (r^2)
    - get single-trial betas
    - get best-fitting HRF
    - get fractional ridge regression regularization level
    - verify reliability and quality control plots manually


## MindEye
1. Transfer derivative outputs from jukebox to Della
    - once GLMsingle has run, move the output (only need TYPED_FITHRF_GLMDENOISE_RR.npz) into the glmsingle directory e.g. /scratch/gpfs/ri4541/MindEyeV2/src/mindeyev2/glmsingle_sub-003_ses-01
    - move the brain mask (..._brain.nii.gz) and the nsdgeneral mask (_nsdgeneral.nii.gz) into the same glmsingle directory as above
    - move the csv file (sub-003_ses-01.csv) to the mindeye directory (home/ri4541/real_time_mindEye2/csv)

2. Run main notebook on Della as a SLURM job
    - update ```/home/ri4541/real_time_mindEye2/run_all_batch.slurm``` with the appropriate variables before submitting
    - ```sbatch /home/ri4541/real_time_mindEye2/run_all_batch.slurm```

3. Sync run with WandB
    - get the wandb sync command from the slurm .err file
    - run that wandb command in Della terminal, example usage: ```wandb sync /home/ri4541/real_time_mindEye2/wandb/offline-run-20250624_162835-sub-005_ses-04_task-B_bs24_MST_rishab_MSTsplit```

4. If desired, run recon_inference, enhanced_recon_inference, and final_evaluations notebooks to get the full list of evaluations

5. Update the [MindEye Evaluations spreadsheet](https://docs.google.com/spreadsheets/d/1-dbmr4ovl2-4-MFNAL1DqLS651KM_ihjDkkUeP1kHXs/edit?usp=sharing) with the new scan