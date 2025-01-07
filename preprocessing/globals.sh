#! /bin/bash

# 1. Edit wherever the dicoms get transferred on the scanner

scanner_dir=/jukebox/dicom/conquest/Prisma-MSTZ400D/NormaL/2024 #Prisma

# 2. Edit where your project directory is
project_dir=/jukebox/norman/rsiyer/rtmindeye

# 3. Edit where your scratch and work directories are located (note: make sure you have setup a work directory on scratch)
scratch_dir=/jukebox/norman/rsiyer/rtmindeye/data/work

data_dir=$project_dir/data
bids_dir=$data_dir/bids
raw_dir=$data_dir/dicom #this is where I want the data from conquest to be copied into my study directory 
defaced_dir=$bids_dir/derivatives/deface #this is where defaced images will end up
scripts_dir=$project_dir/code/preprocessing #directory with my preprocessing scripts, including this one
