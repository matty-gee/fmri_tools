#!/usr/bin/env python

import pandas as pd
import numpy as np
import fmriprep_singularity as fs

# Define your path names
project_dir = '/sc/arion/projects/k23'
bids_root = f'{project_dir}/BIDS_new/'
output_dir = f'{project_dir}/derivatives/'
fs_license = f'{project_dir}/software/license.txt'

# Define your list of subjects (if 'sub-' is at the start of the subject string, it will be removed)
pd_participants = pd.read_csv(f'{project_dir}/code/fmriprep/participants.tsv', delimiter='\t') # add this
participants = { 'participant_id': list(pd_participants['id']) }

# Define the minerva options
image_location = f'{project_dir}/software/fmriprep-20.2.0.simg' # where is the fmriprep-20.2.0.simg file located?
batch_dir = f'{project_dir}/code/batch_dir' # output directory for all batch scripts
minerva_options = {'image_location': image_location,
                   'batch_dir': batch_dir,
                   'project_dir': project_dir}

# Run the fmriprep-docker command through Minerva on the created BIDS directory
fp_singularity = fs.FmriprepSingularityPipeline(participants, bids_root, output_dir, minerva_options, 
                                                queue='Gu', freesurfer=False, cifti_output=False)
fp_singularity.create_singularity_batch()
fp_singularity.run_singularity_batch(participants) # to submit jobs
