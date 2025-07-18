#! /bin/bash

set -e  # Stop on error

module load ants

# ----------------------
# Help Message
# ----------------------
usage() {
    echo ""
    echo "Usage: $0 <subj> <session_label> <task_name> <mask_name> [--ref_session=<session>] [--multisession] [--dry-run] [--data-folder=<folder>]"
    echo ""
    echo "Description: Converts the NSDGeneral mask to functional space, supporting single and multi-session data."
    echo ""
    echo "Positional Arguments:"
    echo "  subj             Subject ID (e.g., '005')"
    echo "  session_label    Session label (e.g., 'ses-01' for single-session or 'ses-01-02' for multi-session)"
    echo "  task_name        Task name (e.g., 'C')"
    echo "  mask_name        Mask name (e.g., '_task-C')"
    echo ""
    echo "Optional Arguments:"
    echo "  --ref_session=<session>  Specify a reference session for alignment (default: first session in label)"
    echo "  --multisession           Enable multi-session processing using fmriprep's anat folder"
    echo "  --dry-run                Perform a dry run (prints actions instead of applying transforms)"
    echo "  --data-folder=<folder>   Specify a custom data folder (default: data_sub-xxx)"
    echo "  --help                   Show this help message and exit"
    echo ""
    echo "Example Usages:"
    echo "  Single-session:       $0 005 ses-01 C _task-C"
    echo "  Multi-session:        $0 005 ses-01-02 C _task-C --multisession"
    echo "  Default data folder:  $0 005 ses-01 C _task-C"
    echo "  Custom data folder:   $0 005 ses-01 C _task-C --data-folder=my_custom_folder"
    echo "  Dry run:              $0 005 ses-01 C _task-C --dry-run"
    echo ""
    exit 0
}

# Check if --help is provided
for arg in "$@"; do
    if [[ "$arg" == "--help" ]]; then
        usage
    fi
done

# ----------------------
# Argument Parsing
# ----------------------
# Required positional arguments:
subj=$1
session_label=$2  # Now takes "ses-01-02" for multi-session cases
task_name=$3
mask_name=$4

# Optional keyword arguments:
ref_session=""   # Default: unset (will be determined dynamically)
multisession=0   # Default: single-session
dry_run=0        # Default: false
data_folder=""   # Default: empty, will be set dynamically

# Parse optional keyword arguments
for arg in "$@"; do
    case $arg in
        --ref_session=*)
            ref_session="${arg#*=}"
            ;;
        --multisession)
            multisession=1
            ;;
        --dry-run)
            dry_run=1
            ;;
        --data-folder=*)
            data_folder="${arg#*=}"
            ;;
        *)
            ;;
    esac
done

# ----------------------
# Directory Setup
# ----------------------
SUBJ_DIR=sub-$subj
STUDY_DIR=/jukebox/norman/rsiyer/rtmindeye
DATA_DIR=$STUDY_DIR/${data_folder:-data_$SUBJ_DIR}  # Use provided folder or default
BIDS_DIR=$DATA_DIR/bids
SCRIPT_DIR=$STUDY_DIR/code/analysis
DERIV_DIR=$BIDS_DIR/derivatives
MASK_DIR=$DERIV_DIR/masks/$SUBJ_DIR

roi_img_path=$SCRIPT_DIR/nsdgeneral_to_MNI.nii.gz

# ----------------------
# Set MNI-to-T1 Transform Path
# ----------------------
if [[ "$session_label" =~ ^ses-[0-9]+$ ]]; then
    # Single-session case
    ref_session=${ref_session:-$session_label}
    mni_to_T1_xfm=$DERIV_DIR/fmriprep/$SUBJ_DIR/${ref_session}/anat/${SUBJ_DIR}_${ref_session}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
    first_session=$session_label
else
    # Multi-session case
    multisession=1
    mni_to_T1_xfm=$DERIV_DIR/fmriprep/$SUBJ_DIR/anat/${SUBJ_DIR}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5

    # Extract first 6 characters (e.g., ses-01 from ses-01-02)
    first_session=${session_label:0:6}
fi

subject_template_path=$DERIV_DIR/fmriprep/$SUBJ_DIR/${first_session}/func/${SUBJ_DIR}_${first_session}_task-${task_name}_run-01_space-T1w_boldref.nii.gz

# ----------------------
# Set Output Mask Filename
# ----------------------
roi_in_epi_space=$MASK_DIR/${SUBJ_DIR}_${session_label}${mask_name}_nsdgeneral.nii.gz

files=("$roi_img_path" "$subject_template_path" "$mni_to_T1_xfm")
for file in "${files[@]}"; do
    if [[ ! -e "$file" ]]; then
        echo "Error: File not found - $file"
        exit 1
    fi
done

# ----------------------
# Apply Transformation
# ----------------------
if [[ $dry_run -eq 1 ]]; then
    echo "Dry Run Mode: The following files will be used:"
    
    for file in "${files[@]}"; do
        if [[ -e "$file" ]]; then
            echo "  [✔] $file (exists)"
        else
            echo "  [✘] $file (MISSING)"
        fi
    done

    echo -e "\nDry Run: The following command would be executed:"
    echo "antsApplyTransforms -i $roi_img_path -r $subject_template_path -t [$mni_to_T1_xfm,0] -n NearestNeighbor -o $roi_in_epi_space -d 3 --float 0 -v 1"
else
    antsApplyTransforms \
        -i $roi_img_path \
        -r $subject_template_path \
        -t [$mni_to_T1_xfm,0] \
        -n NearestNeighbor -o $roi_in_epi_space \
        -d 3 --float 0 -v 1

    echo "NSDGeneral mask saved to: $roi_in_epi_space"
fi
