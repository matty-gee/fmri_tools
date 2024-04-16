#!/bin/bash

#--------------------------------------------------------------------------------------------------
# Run dcm2bids conversion 

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
# - $ sh bids_conversion.sh

## SHOULD REFLECT OWN PATHS
code_dir=/Users/matty_gee/Dropbox/Projects/fmri_tools/preprocessing/dcm2bids/lib 
dcms_dir=/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/Scans/DICOMS
bids_dir=/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/Scans/BIDS

rm -rf $bids_dir/tmp_dcm2bids

# find the number of folders in the dcms_dir
subject_folders=()
for item in $dcms_dir/*; do
  if [ -d "$item" ]; then
    subject_folders+=("$item")
  fi
done

num_subject=`ls -d $dcms_dir/*/ | wc -l`
echo $num_subject " subjects to convert"
 
#--------------------------------------------------------------------------------------------------
# CONVERT TO NIFTI WITH BIDS FORMATTING
#--------------------------------------------------------------------------------------------------

# loop over each folder in subject_folders
for i in ${!subject_folders[@]} ; do

    subject_id=${subject_folders[$i]##*/}
    # mv ${subject_folders[$i]} $dcms_dir/${subject_id}

    #------------- CONVERT TO NIFTI AND PUT IN BIDS -------------#
    
    echo "${i} BIDS conversion for subject" ${subject_id}

    dcm2bids \
      -d $dcms_dir/${subject_id}/dicoms \
      -p ${subject_id} \
      -c $code_dir/dcm2bids_config.json \
      -o $bids_dir \
      #--clobber --forceDcm2niix

    # #------------- MODIFY JSON FILES -------------#
    # # add task name to .json file. eg: "TaskName": "${task}"
    # # add IntendedFor for fmap

    # echo ${subject_id}
    # $code_dir/modify_json.py \
    #   --subject ${subject_id} \
    #   --bids_path $bids_dir

done


#--------------------------------------------------------------------------------------------------
# UNUSED
#--------------------------------------------------------------------------------------------------
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
