#!/bin/bash

#--------------------------------------------------------------------------------------------------
# Download dicoms from XNAT

# requires: 
# - environment.yml : load in python environment w/ req'd python modules: dcm2niix, dcm2bids & xnat
# -- create conda enviornment: $ conda env create -f environment.yml
# -- activate conda environment: $ conda activate dcm2bids
# - lib/dcm2bids_config.json :
# - lib/modify_json.py : 
# - lib/xnat_credentials.txt :
# - lib/xnat_download.py :

# to run from terminal: 
# - $ conda activate dcm2bids
# - $ sh download_dicoms.sh

## UPDATE THESE TO REFLECT OWN PATHS
code_dir=/Users/matty_gee/Dropbox/Projects/preprocessing/dcm2bids/lib 
dcms_dir=/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/Scans/DICOM

#--------------------------------------------------------------------------------------------------
# STEP 1: DOWNLOAD XNAT DICOM DATA
#--------------------------------------------------------------------------------------------------

range_subject="83-84" # specify the index (from 1) to download dcms for - need a better way....
num_subject=1 # how many subjects to download

echo "Start downloading XNAT data for" $num_subject "participant(s)"
$code_dir/xnat_download.py \
  --range_subject $range_subject \
  --username `sed -n 1p $code_dir/xnat_credentials.txt` \
  --password `sed -n 2p $code_dir/xnat_credentials.txt` \
  --download_dir $dcms_dir