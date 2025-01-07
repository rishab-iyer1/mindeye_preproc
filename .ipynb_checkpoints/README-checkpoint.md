# Data Processing Pipeline
## Pre-processing
Start with raw DICOMs, repeat following steps for each scanning session of a subject

---

1. ```preprocessing/step1_preproc.sh```
    - use heudiconv to convert DICOMs to BIDS format
    - use PyDeface to deface images
      
2. ```preprocessing/step2_preproc.sh```
    - delete scouts and duplicate runs
    - modify fieldmap files to be BIDS-compatible
    - current issues
        - IntendedFor field in fmap folder gets added repeatedly if it already exists causing JSON errors
        - .tsv file doesn't update with deleted files
      
3. [BIDS Validator](https://bids-standard.github.io/bids-validator/)
    - verify manually that the data is BIDS-compatible
      
4. ```preprocessing/run_mriqc.sh```
    - quality control checks, verify outputs manually
      
5. ```preprocessing/run_mriqc_group.sh```
    - quality control checks, verify outputs manually

---

1. ```preprocessing/run_fmriprep.sh```
    - runs fMRIPrep for all sessions for an individual participant

## Analysis
1. ```analysis/nsdgeneral_to_epi.sh```
    - create an NSDgeneral mask for a subject containing primarily visual cortex voxels
      
2. ```analysis/GLMsingle.ipynb```
    - get model accuracy (r^2)
    - get single-trial betas
    - get best-fitting HRF
    - get fractional ridge regression regularization level
    - verify reliability and quality control plots manually