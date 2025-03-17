#!/usr/bin/env python3

import time
from tqdm import tqdm
import subprocess
from datetime import datetime

def test_tqdm_logging():
    # Create a log file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"tqdm_test_{timestamp}.log"
    
    print(f"Running tqdm test...")
    print(f"Logs will be saved to: {log_file}")
    
    start_time = time.time()
    
    with open(log_file, 'w') as f:
        # Write header
        f.write(f"TQDM Test Log\n")
        f.write(f"============\n")
        f.write(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Test 1: Simple progress bar
        f.write("Test 1: Simple progress bar\n")
        f.write("-" * 20 + "\n")
        for i in tqdm(range(5), desc="Simple progress", file=f):
            time.sleep(1)
            f.flush()
        
        # Test 2: Nested progress bars
        f.write("\nTest 2: Nested progress bars\n")
        f.write("-" * 20 + "\n")
        for i in tqdm(range(3), desc="Outer loop", file=f):
            for j in tqdm(range(2), desc=f"Inner loop {i}", file=f, leave=False):
                time.sleep(0.5)
                f.flush()
        
        # Test 3: Progress bar with additional output
        f.write("\nTest 3: Progress bar with additional output\n")
        f.write("-" * 20 + "\n")
        pbar = tqdm(range(3), desc="Progress with output", file=f)
        for i in pbar:
            time.sleep(1)
            f.write(f"Additional output at step {i}\n")
            f.flush()
            pbar.set_postfix({"step": i})
    
    end_time = time.time()
    duration = end_time - start_time
    
    # Write execution summary
    with open(log_file, 'a') as f:
        f.write(f"\nExecution Summary\n")
        f.write(f"================\n")
        f.write(f"End Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Duration: {duration:.2f} seconds\n")
    
    print(f"\nTest complete. Logs saved to: {log_file}")
    print(f"Total Duration: {duration:.2f} seconds")

if __name__ == "__main__":
    test_tqdm_logging() 