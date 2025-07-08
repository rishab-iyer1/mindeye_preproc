#!/usr/bin/env python3

import os
import argparse
import subprocess
import time
from datetime import datetime

def run_glmsingle(data_dir, glmsingle_dir):
    """
    Run GLMsingle analysis with proper directory and logging setup
    
    """
    # Set up directory and logging
    # glmsingle_dir, log_file = setup_glmsingle_dir(sub, ses_list, task_label=task_label)
    assert os.path.exists(f"/usr/people/ri4541/rtmindeye/{data_dir}"), f'data directory does not exist:\n{data_dir}'
    for term in ("glmsingle_", "_ses-", "_task-"):
        assert term in str(glmsingle_dir)
    output_dir = f"/usr/people/ri4541/rtmindeye/{data_dir}/bids/derivatives/{glmsingle_dir}"
    if not os.path.exists(output_dir):
        ask = input(f'The specified glmsingle path does not exist:\n{output_dir}\nWould you like to create this path? (y/n): ')
        if ask in ('Y', 'y'):
            os.makedirs(output_dir)
        elif ask in ('N', 'n'):
            raise FileNotFoundError
        else:
            raise ValueError('Invalid response. Please respond with "y" or "n".')
    # err
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f"execution_{timestamp}.log")

    # Convert notebook to script
    notebook_path = '/usr/people/ri4541/rtmindeye/code/analysis/GLMsingle.ipynb'
    script_path = os.path.join(output_dir, 'GLMsingle.py')
    
    print(f"Converting notebook to script...")
    subprocess.run([
        'jupyter', 'nbconvert',
        '--to', 'script',
        notebook_path,
        '--output-dir', output_dir
    ], check=True)
    
    # Run GLMsingle with logging
    print(f"Running GLMsingle analysis...")
    print(f"Logs will be saved to: {log_file}")
    
    start_time = time.time()
    
    with open(log_file, 'w') as f:
        # Write header with execution info
        f.write(f"GLMsingle Execution Log\n")
        f.write(f"=====================\n")
        f.write(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Run the script and capture output
        process = subprocess.Popen(
                    ['ipython', script_path],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    bufsize=1  # Line buffering
                )
        # Write output to both console and log file
        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line:
                print(line, end='', flush=True)
                f.write(line)
                f.flush()  # Ensure immediate writing to log file
    
    end_time = time.time()
    duration = end_time - start_time
    hours = int(duration // 3600)
    minutes = int((duration % 3600) // 60)
    seconds = int(duration % 60)
    
    # Write execution summary
    with open(log_file, 'a') as f:
        f.write(f"\nExecution Summary\n")
        f.write(f"================\n")
        f.write(f"End Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Duration: {hours:02d}:{minutes:02d}:{seconds:02d}\n")
    
    print(f"\nAnalysis complete. Logs saved to: {log_file}")
    print(f"Total Duration: {hours:02d}:{minutes:02d}:{seconds:02d}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GLMsingle with specified config")
    parser.add_argument("data_dir", help="Data directory (i.e. data_sub-005_ses-04)")
    parser.add_argument("glmsingle_dir", help="GLMsingle directory (i.e. glmsingle_ses-04_task-B_resampled_2_0mm_trilinear)")
    parser.add_argument("sub", type=str, help="Subject ID (e.g., 'sub-001')")
    parser.add_argument("session", type=str, help="Session ID (e.g., 'ses-01' (single-session), 'all' (multi-session))")
    parser.add_argument("func_task_name", type=str, help="Functional task name (e.g., 'study')")
    parser.add_argument("--resample_voxel_size", action='store_true', default=False, help="Resample voxel size flag (True/False)")
    parser.add_argument("--run_resample_voxel", action='store_true', default=True, help="Whether to perform resampling (default: True). If False, assumes resampled data has been previously computed and loads it in. If True but resample_voxel_size is False, has no effect.")
    parser.add_argument("--resampled_vox_size", type=float, help="Voxel size (mm isotropic) to resample to (e.g., 2.5)")
    parser.add_argument("--resample_method", type=str, help="Resampling method (e.g., 'trilinear')")
    parser.add_argument("--ses_list", type=str, help="Comma-separated list of session IDs (e.g. ses-01,ses-02).")
    parser.add_argument("--design_ses_list", type=str, help="list of design matrix session IDs (e.g. ['ses-01', 'ses-02']). Use only if session='all' and if the design matrix session ID doesn't match the scan's session ID.")
    parser.add_argument("--ref_session", type=str, help="Reference session to draw anatomical from (e.g., 'ses-xx'). Use this only if the current session being analyzed does not have a T1 of its own. Make sure the current session and the provided reference session were fMRIPrepped together.")
    parser.add_argument('--dry_run', action='store_true', help="Run script without executing GLMsingle.")

    args = parser.parse_args()
    
    # argument validation
    # TODO 

    # set all args to be env variables; to be read in when running the main script
    for arg, value in vars(args).items():
        if args.dry_run:
            print(f"os.environ[{arg.upper()}] -> {str(value)}")
        else:
            os.environ[arg.upper()] = str(value) if value is not None else ""
    
    if not args.dry_run:
        run_glmsingle(args.data_dir, args.glmsingle_dir)