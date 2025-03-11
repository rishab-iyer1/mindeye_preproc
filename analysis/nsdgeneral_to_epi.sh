#! /bin/bash

set -e #stop immediately when error occurs
# example usage: ./nsdgeneral_to_epi.sh 001 02 C _task-C
# where C is the task_name and _task-C is the mask_name (task_name and mask_name are optional)
# multi-session: ./nsdgeneral_to_epi.sh 005 01 C _task-C all --multisession


module load ants

# required positional arguments
subj=$1
session=$2
task_name=$3
mask_name=$4

# optional keyword arguments
ref_session=$session  # Default to session
multisession=0  # Default to single-session

# Parse optional keyword arguments
for arg in "$@"; do
    case $arg in
        --ref_session=*)
            ref_session="${arg#*=}"
            ;;
        --multisession)
            multisession=1
            ;;
        *)
            ;;
    esac
done

SUBJ_DIR=sub-$subj
STUDY_DIR=/jukebox/norman/rsiyer/rtmindeye
DATA_DIR=$STUDY_DIR/data_$SUBJ_DIR
BIDS_DIR=$DATA_DIR/bids
SCRIPT_DIR=$STUDY_DIR/code/analysis
DERIV_DIR=$BIDS_DIR/derivatives
FREESURFER_DIR=$DERIV_DIR/sourcedata/freesurfer/$SUBJ_DIR/mri
MASK_DIR=$DERIV_DIR/masks/$SUBJ_DIR

roi_img_path=$SCRIPT_DIR/nsdgeneral_to_MNI.nii.gz

if [ "$multisession" = 1 ]; then
    mni_to_T1_xfm=$DERIV_DIR/fmriprep/$SUBJ_DIR/anat/${SUBJ_DIR}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
else
    mni_to_T1_xfm=$DERIV_DIR/fmriprep/$SUBJ_DIR/ses-${ref_session}/anat/${SUBJ_DIR}_ses-${ref_session}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
fi
subject_template_path=$DERIV_DIR/fmriprep/$SUBJ_DIR/ses-${session}/func/${SUBJ_DIR}_ses-${session}_task-${task_name}_run-01_space-T1w_boldref.nii.gz
roi_in_epi_space=$MASK_DIR/${SUBJ_DIR}_ses-${session}${mask_name}_nsdgeneral.nii.gz

FILES=("$roi_img_path" "$mni_to_T1_xfm" "$subject_template_path" "$roi_in_epi_space")

# Loop through the list
for FILE in "${FILES[@]}"; do
    if [ -e "$FILE" ]; then
        echo "$FILE exists."
    else
        echo "$FILE does not exist."
    fi
done

# echo $roi_in_epi_space

# antsApplyTransforms \
#         -i $roi_img_path \
#         -r $subject_template_path \
#         -t [$mni_to_T1_xfm,0]  \
#         -n NearestNeighbor -o $roi_in_epi_space \
#         -d 3 --float 0 \
#         -v 1