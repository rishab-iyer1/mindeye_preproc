"""
Utility functions to prepare data for GLMsingle
"""

import os
from functools import wraps

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.maskers import NiftiMasker
from nilearn.plotting import plot_roi

def get_session_label(ses_list):
    """
    Create session label from list of session strings
    
    Args:
        ses_list: List of session strings e.g. ['ses-01', 'ses-02']
    
    Returns:
        session_label: e.g. 'ses-01-02' or 'ses-01+03'
    
    Examples:
        ['ses-01', 'ses-02', 'ses-03'] -> 'ses-01-03'
        ['ses-01', 'ses-03'] -> 'ses-01+03'
    """
    assert isinstance(ses_list, list), "ses_list must be a list"
    assert all(isinstance(ses, str) for ses in ses_list), "All elements in ses_list must be strings"
    assert all(ses.startswith('ses-') for ses in ses_list), "All elements in ses_list must start with 'ses-'"
    assert len(ses_list) > 0, "ses_list must not be empty"

    # Extract numbers from session strings
    sessions = [int(ses.split('-')[1]) for ses in ses_list]
    
    if len(sessions) == 1:
        return ses_list[0]
        
    # Sort sessions
    sessions = sorted(sessions)
    
    # Find consecutive ranges
    ranges = []
    range_start = sessions[0]
    prev = sessions[0]
    
    for curr in sessions[1:] + [None]:
        if curr != prev + 1:
            range_end = prev
            if range_start == range_end:
                ranges.append(f"{range_start:02d}")
            else:
                ranges.append(f"{range_start:02d}-{range_end:02d}")
            range_start = curr
        prev = curr
    
    return "ses-" + "+".join(ranges)


def process_design(filename):
    """
    Process design CSV file and extract relevant trial information
    
    Args:
        filename: Path to CSV file containing design information
    
    Returns:
        data: Pandas DataFrame containing all design information
        starts: Array of trial start times
        images: Array of image names
        is_new_run: Array of boolean flags indicating new runs
    """
    data = pd.read_csv(filename)
    data = data.dropna(subset=['current_image'])  # in these csv files, there are blanks between runs
    starts = data['trial.started'].values[10:]
    images = data['current_image'].values[10:]
    is_new_run = data['is_new_run'].values[10:]
    return data, starts, images, is_new_run


def create_design_matrix(images, starts, is_new_run, unique_images, n_runs, n_trs, len_unique_images):
    """
    Create a design matrix for GLMsingle analysis
    
    Args:
        images: List of image filenames
        starts: List of trial start times
        is_new_run: List of boolean flags indicating if trial is last in run
        unique_images: List of unique image filenames
        n_runs: Number of runs
        n_trs: Number of TRs per run
        len_unique_images: Length of unique_images list
    
    Returns:
        design: List of design matrices, one per run
    """
    design = [np.zeros((n_trs, len_unique_images)) for _ in range(n_runs)]
    starting_time = starts[0]
    cur_run = 0
    first_trial_of_new_run = False
    
    for i in range(len(images)):
        if first_trial_of_new_run: # is the first trial of a new run?
            starting_time = starts[i]
            cur_run += 1
            first_trial_of_new_run = False
            
        if is_new_run[i] == 1: # is this the last trial of the run?
            first_trial_of_new_run = True
            if n_runs == 1:
                break
        
        if images[i] == "blank.jpg":
            continue
    
        image_idx = np.where(str(images[i])==unique_images)[0].item()
        timepoint = int(np.round(starts[i] - starting_time))
        
        design[cur_run][timepoint, image_idx] = 1
            
    return design


def mask_info(avg_mask, t1, plot=False):
    """
    Print information about the brain mask
    
    Args:
        avg_mask: nibabel.Nifti1Image object of the average brain mask
        t1: nibabel.Nifti1Image object of the T1-weighted image
        plot: boolean flag to plot the brain mask
    
    Returns:
        dimsize: Tuple of voxel dimensions
        affine_mat: 4x4 affine transformation matrix
        brain: 3D numpy array of brain mask
        xyz: Tuple of brain dimensions
    """
    if plot:
        plot_roi(avg_mask, bg_img=t1)
    
    # mask info
    dimsize = avg_mask.header.get_zooms()
    affine_mat = avg_mask.affine
    brain = avg_mask.get_fdata()
    xyz = brain.shape  # xyz dimensionality of brain mask and epi data
    
    print('Mask dimensions:', dimsize)
    print('')
    print('Affine:')
    print(affine_mat)
    print('')
    print(f'There are {np.sum(brain)} voxels in the included brain mask\n')

    return dimsize, affine_mat, brain, xyz


def fit_save_bold(in_file, mask, out_file, run):
    """
    Fit and save BOLD data
    
    Args:
        in_file: str, path to the input BOLD file
        mask: nibabel.Nifti1Image object of the brain mask
        out_file: str, path to the output npy file
        run: int, run number
    
    Returns:
        epi_mask_data: numpy array of masked BOLD data
    """
    epi_masker = NiftiMasker(mask_img=mask)
    print('loading data for run', run+1, ':', in_file, '\n')
    
    epi_mask_data = epi_masker.fit_transform(in_file)
    epi_mask_data = epi_mask_data.T  # transpose to make it voxels,time
    
    # save individual run data to npy file
    print(f'\n*** Saving data to {out_file} ***\n')
    np.save(out_file, epi_mask_data)
    return epi_mask_data


def compare_mask_epi_dims(epi_file, mask):
    """
    Prints info about epi data and compares to input mask

    Args:
        epi_file: str, file path to epi data
        mask: nifti object containing the mask to compare
    
    Returns:
        epi_dimsize: Tuple of EPI voxel dimensions
        epi_affine: 4x4 affine transformation matrix
    """
    # get some info about epi data by loading study run 1
    epi_data = nib.load(epi_file)
    run1 = epi_data.get_fdata()
    
    print('checking epi data for run 1:', epi_file)
    print('')
    epi_dimsize = epi_data.header.get_zooms()
    epi_affine = epi_data.affine
    print('Dimensions:', epi_dimsize)
    print('Affine:')
    print(epi_data.affine)
    print('')
    
    # get shape of data volume (XYZ) for convenience
    xyz = run1.shape[:3]
    xyzt = run1.shape
    print(xyz)
    print(xyzt)
    
    # double check that brain mask and epi data have same dimensions and affine
    assert mask.header.get_zooms() == epi_dimsize[:3]
    assert mask.affine.all() == epi_data.affine.all()
    
    return epi_dimsize, epi_affine


def log_io(func):
    """
    Decorator that prints loading and saving information for functions
    that process input and output files
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        inp = args[0]
        output = kwargs['output']
        print(f'\n*** Loading data from {inp} ***\n')
        result = func(*args, **kwargs)
        print(f'\n*** Saved resampled data to {output} ***\n')
        return result
    return wrapper


@log_io
def resample(inp, ref, target_size, omat, output=None):
    """
    Resample an image using flirt
    
    Args:
        inp: str, path to the input image
        ref: str, path to the reference image
        target_size: float, target voxel size
        omat: str, path to the output matrix file
        output: str, path to the output image file
    """
    os.system(f"flirt -in {inp} \
                    -ref {ref} \
                    -applyisoxfm {target_size} -nosearch \
                    -omat {omat} \
                    -out {output}")


@log_io
def applyxfm(inp, ref, init, interp, output=None):
    """
    Apply a transformation matrix to an image
    
    Args:
        inp: str, path to the input image
        ref: str, path to the reference image
        init: str, path to the initial transformation matrix
        interp: str, interpolation method
        output: str, path to the output image file
    """
    os.system(f"flirt -in {inp} \
                -ref {ref} \
                -out {output} \
                -applyxfm -init {init} \
                -interp {interp}")


def apply_thresh(inp, thresh, output=None):
    """
    Apply a threshold to an image
    
    Args:
        inp: str, path to the input image
        thresh: float, threshold value
        output: str, path to the output image file
    """
    os.system(f"fslmaths {inp} -thr {thresh} -bin {output}")