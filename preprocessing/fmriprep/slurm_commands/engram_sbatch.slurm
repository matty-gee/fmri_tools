#!/bin/bash
#
#SBATCH --job-name=fmriprep
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=general
# Outputs ----------------------------------
#SBATCH -o /mnt/synapse/projects/SocialSpace/Data/Scans/RT/fmriprep_/code/singularity_log/%x-%A-%a.out
#SBATCH -e /mnt/synapse/projects/SocialSpace/Data/Scans/RT/fmriprep_/code/singularity_log/%x-%A-%a.err
# ------------------------------------------

STUDY=/mnt/synapse/projects/SocialSpace/Data/Scans/RT/fmriprep_
BIDS_DIR="$STUDY/BIDS/"
DERIVS_DIR="$STUDY/derivatives/"
CODE_DIR="$STUDY/code"
BATCH_DIR="$STUDY/batch_dir"
PREP_DIR=/mnt/synapse/projects/SocialSpace/Data/Scans/fmriprep_pipeline
PARTICIPANTS_FILENAME="participants.tsv"

# Prepare some writeable bind-mount points.
TEMPLATEFLOW_HOST_HOME=$CODE_DIR/.cache/templateflow
mkdir -p ${TEMPLATEFLOW_HOST_HOME}

# Designate a templateflow bind-mount point
export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"
SINGULARITY_CMD="singularity run --cleanenv \
				-B $BIDS_DIR:/data \
				-B ${TEMPLATEFLOW_HOST_HOME}:${SINGULARITYENV_TEMPLATEFLOW_HOME} \
				-B $CODE_DIR:/code \
				-B $DERIVS_DIR:/out \
				-B $BATCH_DIR:/batch \
				$PREP_DIR/fmriprep-20.2.0.simg"

# Parse the participants.tsv file and extract one subject ID from the line corresponding to this SLURM task.
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID}))s/sub-(\S*)\>.*/\1/gp" ${CODE_DIR}/${PARTICIPANTS_FILENAME} )
echo $subject

# fMRI prep arguments 
cmd="${SINGULARITY_CMD} /data /out participant \
	--participant-label $subject \
	-w /batch/fmriprep_work \
	-t socialnav \
	--use-aroma \
	--skip_bids_validation \
	--notrack \
	--fs-license-file /code/license.txt \
	--fs-no-reconall"	

# Setup done, run the command
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
eval $cmd
exitcode=$?

# Output results to a table
echo "sub-$subject   ${SLURM_ARRAY_TASK_ID}    $exitcode" \
      >> ${SLURM_JOB_NAME}.${SLURM_ARRAY_JOB_ID}.tsv
echo Finished tasks ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode