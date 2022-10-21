#!/bin/bash

# (1) Download dicoms from XNAT & (2) run dcm2bids conversion 
# By Tien Tong, with edits by Matthew Schafer, 2022

# requires: 
# - environment.yml : load in python environment w/ req'd python modules: dcm2niix, dcm2bids & xnat
# -- create conda enviornment: $ conda env create -f environment.yml
# -- activate conda environment: $ conda activate dcm2bids
# - lib/dcm2bids_config.json :
# - lib/modify_json.py : 
# - lib/xnat_credentials.txt :
# - lib/xnat_download.py :

# to run from terminal: $ sh bids_conversion.sh

## UPDATE THESE TO REFLECT OWN PATHS
code_dir=/Users/matty_gee/Dropbox/Projects/fmri_tools/dcm2bids/lib 
dcms_dir=/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/Scans/DICOM
bids_dir=/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/Scans/BIDS

###############################################################################
####################### STEP 1: DOWNLOAD XNAT DICOM DATA ######################
###############################################################################

# CHANGE THE NEXT TWO LINES:
range_subject="82-83" # specify the index (from 1) to download dcms for
num_subject=4 # how many subjects to convert

rm -rf $bids_dir/tmp_dcm2bids

echo "Start downloading XNAT data for" $num_subject "participant(s)"
$code_dir/xnat_download.py \
  --range_subject $range_subject \
  --username `sed -n 1p $code_dir/xnat_credentials.txt` \
  --password `sed -n 2p $code_dir/xnat_credentials.txt` \
  --download_dir $dcms_dir
 
###############################################################################
################### STEP 2: CONVERT TO NIFTI AND PUT IN BIDS ##################
###############################################################################

new_download=($(ls -td $dcms_dir/* | head -${num_subject}))
for i in ${!new_download[@]} ; do
    subject_id=`echo ${new_download[$i]} | awk -F "/" '{print $NF}'`

    mv ${new_download[$i]} $dcms_dir/${subject_id}

    ################### STEP 2.1a: CONVERT TO NIFTY AND PUT IN BIDS ############
    echo "Start doing BIDS conversion for subject" ${subject_id}
    dcm2bids \
      -d $dcms_dir/${subject_id}/*/scans \
      -p ${subject_id} \
      -c $code_dir/dcm2bids_config.json \
      -o $bids_dir \
      --clobber --forceDcm2niix

    ###############################################################################
    ############################ STEP 3: MODIFY JSON FILES ########################
    ###############################################################################

    # add task name to .json file. eg: "TaskName": "${task}"
    # add IntendedFor for fmap

    echo ${subject_id}
    $code_dir/modify_json.py \
      --subject ${subject_id} \
      --bids_path $bids_dir

done

########################################################################################################
##### UNUSED
########################################################################################################
#     # ######################## STEP 2.1b: RENAME FILES ###########################
#     # # run and echo were flipped using dcm2bids, have to switch them back
#     # subject_func=$bids_dir/sub-${subject_id}/func

#     # for task in rest whyhow ; do
#     #   task_list=($(ls $subject_func/*${task}*bold.nii.gz))
#     #   run_num=$((${#task_list[@]} / 3))
#     #   if [ $run_num -gt 1 ] ; then
#     #   for run in `seq -s " " -f "%02g" $run_num` ; do
#     #     for echo in 01 02 03 ; do
#     #       for type in bold sbref ; do
#     #         for extension in json nii.gz ; do
#     #           old=$subject_func/sub-${subject_id}_task-${task}_echo-${echo}_run-${run}_${type}.${extension}
#     #           new=$subject_func/sub-${subject_id}_task-${task}_run-${run}_echo-${echo}_${type}.${extension}
#     #           mv $old $new
#     #         done
#     #       done
#     #     done
#     #   done
#     #   fi
#     # done

#     # echo "BIDS conversion completed for subject" ${subject_id}

#     # ################### STEP 2.2: RUN OPTIMAL TE COMBINATION #####################
#     # for task in rest whyhow ; do
#     #   task_list=($(ls $subject_func/*${task}*bold.nii.gz))
#     #   run_num=$((${#task_list[@]} / 3))
#     #   echo "Start the tedana optimal combination for" $run_num "run(s) of the" ${task} "task of subject" ${subject_id}

#     #   if [ $run_num -gt 1 ] ; then
#     #     for run in `seq -s " " -f "%02g" $run_num` ; do
#     #       input=(`ls $subject_func/*${task}*run-${run}*bold.nii.gz`)
#     #       t2smap \
#     #         -d ${input[@]} \
#     #         -e 10.80 28.68 46.56 \
#     #         --out-dir $subject_func
#     #       rm $subject_func/*S0map* $subject_func/*T2starmap*
#     #       mv $subject_func/desc-optcom_bold.nii.gz $subject_func/sub-${subject_id}_task-${task}_run-${run}_bold.nii.gz
#     #       cp $subject_func/sub-${subject_id}_task-${task}_run-${run}*echo-01*_bold.json $subject_func/sub-${subject_id}_task-${task}_run-${run}_bold.json
#     #     done
#     #   else
#     #     input=(`ls $subject_func/*${task}*bold.nii.gz`)
#     #     t2smap \
#     #       -d ${input[@]} \
#     #       -e 10.80 28.68 46.56 \
#     #       --out-dir $subject_func
#     #     rm $subject_func/*S0map* $subject_func/*T2starmap*
#     #     mv $subject_func/desc-optcom_bold.nii.gz $subject_func/sub-${subject_id}_task-${task}_bold.nii.gz
#     #     cp $subject_func/sub-${subject_id}_task-${task}*echo-01*_bold.json $subject_func/sub-${subject_id}_task-${task}_bold.json

#     #   fi
#       # need to delete echo time of the optimally combined bold
#     # done
