#!/usr/bin/env python3

import os
import sys
import subprocess
import time
from datetime import datetime
# from utils import setup_glmsingle_dir

def run_glmsingle(sub, session, ses_list, data_dir, glmsingle_dir):
    """
    Run GLMsingle analysis with proper directory and logging setup
    
    Args:
        sub: str, subject ID (e.g. 'sub-001')
        session: str, session ID (e.g. 'ses-02' or 'all')
        ses_list: list, list of session IDs
    """
    # Set up directory and logging
    # glmsingle_dir, log_file = setup_glmsingle_dir(sub, ses_list, task_label=task_label)
    output_dir = f"/usr/people/ri4541/rtmindeye/{data_dir}/bids/derivatives/{glmsingle_dir}"

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f"execution_{timestamp}.log")

    # Convert notebook to script
    notebook_path = os.path.join(os.path.dirname(__file__), 'GLMsingle.ipynb')
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
        f.write(f"Subject: {sub}\n")
        f.write(f"Session: {session}\n")
        f.write(f"Session List: {ses_list}\n")
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
    if len(sys.argv) < 3:
        print("Usage: python run_glmsingle.py <subject> <session> <data_dir> <glmsingle_dir>")
        print("Example: python ~/rtmindeye/code/analysis/run_glmsingle.py sub-001 ses-05 data_sub-001_ses-05 glmsingle_task-C_resampled_2_5mm_sinc")
        sys.exit(1)
    
    sub = sys.argv[1]
    session = sys.argv[2]
    data_dir = sys.argv[3]
    glmsingle_dir = sys.argv[4]
    
    # Set up session list based on input
    if session == "all":
        ses_list = ["ses-02", "ses-03"]
    else:
        ses_list = [session]
    # print(ses_list)
    run_glmsingle(sub, session, ses_list, data_dir=data_dir, glmsingle_dir=glmsingle_dir) 