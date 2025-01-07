#! /bin/bash

set -e #stop immediately when error occurs
# example usage: ./nsdgeneral_to_epi.sh 001 02

module load ants

subj=$1
session=$2

SUBJ_DIR=sub-$subj
STUDY_DIR=/jukebox/norman/rsiyer/rtmindeye
DATA_DIR=$STUDY_DIR/data
BIDS_DIR=$DATA_DIR/bids
SCRIPT_DIR=$STUDY_DIR/code/analysis
DERIV_DIR=$BIDS_DIR/derivatives
FREESURFER_DIR=$DERIV_DIR/sourcedata/freesurfer/$SUBJ_DIR/mri
MASK_DIR=$DERIV_DIR/masks/$SUBJ_DIR

# T1=$DERIV_DIR/fmriprep/$SUBJ_DIR/anat/${SUBJ_DIR}_desc-preproc_T1w.nii.gz

roi_img_path=$SCRIPT_DIR/nsdgeneral_to_MNI.nii.gz

# Warp it to subject space
mni_to_T1_xfm=$DERIV_DIR/fmriprep/$SUBJ_DIR/ses-${session}/anat/${SUBJ_DIR}_ses-${session}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
# mni_to_T1_xfm=$DERIV_DIR/fmriprep/$SUBJ_DIR/anat/${SUBJ_DIR}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
subject_template_path=$DERIV_DIR/fmriprep/$SUBJ_DIR/ses-${session}/func/${SUBJ_DIR}_ses-${session}_task-study_run-01_space-T1w_boldref.nii.gz
roi_in_epi_space=$MASK_DIR/${SUBJ_DIR}_nsdgeneral.nii.gz

antsApplyTransforms \
        -i $roi_img_path \
        -r $subject_template_path \
        -t [$mni_to_T1_xfm,0]  \
        -n NearestNeighbor -o $roi_in_epi_space \
        -d 3 --float 0 \
        -v 1