#!/usr/bin/env python
# Author: Kaustubh Kulkarni
# Date: Feb 20, 2020

# For fmriprep-setup
import os, glob, shutil, sys, subprocess, logging, datetime, time, json
import numpy as np
from copy import deepcopy

class FmriprepSingularityPipeline(object):
    """ 
    This class prepares batch scripts and runs fmriprep through a Singularity image.
    Designed for use on Minerva at Mount Sinai.
  
    Note that the following methods must be run in order:
    1. create_singularity_batch()
    2. run_singularity_batch()

    Alternatively, you can run only create_singularity_batch() and submit generated scripts manually. 
  
    """

    def __init__(self, participants, bids_root, output, minerva_options, queue='Gu', freesurfer=False, cifti_output=False):
        """ 
        Constructs the necessary attributes for the FmriprepSingularityPipeline instance.   
      
        Parameters: 
        participants (list): Dictionary of participant IDs and related properties
        bids_root (str): Path to generated BIDS directory
        output (str): Path to store fmriprep output
        minerva_options (dict): Dictionary of HPC-specific options (see below)
        freesurfer (bool): Flag to specify Freesurfer surface estimation
      
        minerva_options dictionary should contain:
            1. image_location: Path to directory that contains the fmriprep singularity image and Freesurfer license file
            2. batch_dir: Path to a directory to store generated batch scripts
            3. project_dir: Root level project directory (parent to bids_root)

        Returns: 
        obj: FmriprepSingularityPipeline instance 
      
        """

        # Define class variables
        self.subs = participants['participant_id']
        self.bids_root = bids_root
        self.output = output
        self.freesurfer = freesurfer
        self.minerva_options = minerva_options
        self.batch_dir = minerva_options['batch_dir']
        self.cifti_output = cifti_output
        self.queue = queue

        if cifti_output and not freesurfer:
            logging.error('Freesurfer must be on to have cifti-output!')
            raise OSError('Freesurfer must be on to have cifti-output!')

    def create_singularity_batch(self):
        """ 
        Creates the subject batch scripts for running fmriprep with Singularity.
        To run in parallel, subjects are run individually and submitted as separate jobs on the cluster.   
      
        """

        logging.info('Setting up fmriprep command through Singularity for Minerva')
        
        # different node options: to make sure am being fair w/ resources
        if self.queue == 'Gu':
            project_acc = 'acc_guLab'
            queue       = 'private'
            self.sleep_time = 3600 * 5 # 3,600s = 1h
        else:
            project_acc = 'acc_k23'
            queue       = self.queue
            self.sleep_time = 60

        # Check if the singularity image exists in the image location
        if not os.path.isfile(f'{self.minerva_options["image_location"]}'):
            logging.error('fmriprep image does not exist in the given location!')
            raise OSError('fmriprep image does not exist in the given directory!')

        # Create the specified batch directory folder if it doesn't exist
        logging.info('Setting up batch directory for subject scripts')
        if not os.path.isdir(self.batch_dir):
            os.makedirs(self.batch_dir)
            logging.info(f'Creating batch dir at {self.batch_dir}')
        else:
            logging.info(f'{self.batch_dir} exists!')
        
        # Create the batch_output folder within the batch directory folder
        # This will hold the outputs from each of the batch scripts
        if not os.path.isdir(f'{self.batch_dir}/batch_output'):
            os.makedirs(f'{self.batch_dir}/batch_output')

        # Loop over all subjects
        for sub in self.subs:

            # Strip the 'sub-' prefix from the subject name string, if it's there
            if sub[:4] == 'sub-': sub = sub[4:]

            # Create the subject specific batch script
            sub_batch_script = f'{self.batch_dir}/sub-{sub}.sh'
            with open(sub_batch_script, 'w') as f:
                lines = [ # These are the BSUB cookies
                    f'#!/bin/bash\n\n',
                    f'#BSUB -J fmriprep_sub-{sub}\n',
                    f'#BSUB -P {project_acc}\n', 
                    f'#BSUB -q {queue}\n', 
                    f'#BSUB -n 3\n',
                    f'#BSUB -W 12:00\n',
                    f'#BSUB -R rusage[mem=16000]\n',
                    f'#BSUB -o {self.batch_dir}/batch_output/nodejob-fmriprep-sub-{sub}.out\n',
                    f'#BSUB -L /bin/bash\n\n',
                    f'ml singularity/3.6.4\n\n', # Module load singularity
                    f'cd {self.minerva_options["project_dir"]}\n\n', # Enter the directory that contains the fmriprep.20.0.1.simg
                    # f'export SINGULARITYENV_TEMPLATEFLOW_HOME=$HOME/.cache/templateflow\n\n' # for template flow on minerva
                ]
                f.writelines(lines)

                # Create the command
                # -- can also ignore slicetiming if necessary...
                command = f"singularity run -B $HOME:/home --home /home \
                            -B {os.path.dirname(self.minerva_options['image_location'])}:/software \
                            --cleanenv \
                            {self.minerva_options['image_location']} \
                            {self.bids_root} \
                            {self.output} \
                            participant \
                            --output-spaces MNI152NLin2009cAsym:res-2 \
                            --participant-label {sub} \
                            -w /software/fmriprep_work \
                            -t socialnav \
                            --use-aroma \
                            --ignore fieldmaps \
                            --skip_bids_validation \
                            --notrack \
                            --fs-license-file /software/license.txt"
                command = " ".join(command.split())
                # Ignore freesurfer if specified
                if not self.freesurfer:
                   command = " ".join([command, '--fs-no-reconall'])
                # Add cifti-output if specified
                if self.cifti_output:
                    command = " ".join([command, '--cifti-output'])
                # Output command to batch script
                f.write(command)

        # Include all variables in the 'minerva_option' dictionary
        self.minerva_options['subs']       = self.subs
        self.minerva_options['bids_root']  = self.bids_root
        self.minerva_options['output']     = self.output
        self.minerva_options['freesurfer'] = self.freesurfer

        # Save all parameters within the batch directory as well
        with open(f'{self.batch_dir}/minerva_options.json', 'w') as f:
            json.dump(self.minerva_options, f) 

    def run_singularity_batch(self, overwrite=False):
        """ 
        Submits generated subject batch scripts to the HPC. 
        Parameters: 
        - subs (list): A list of subject ID strings. May be a subset of subjects in the BIDS directory.
        """
        logging.info('Submitting singularity batch scripts to the queue')
        counter = 1
        for sub in self.subs:
            # Submit job to scheduler
            if sub.startswith('sub-'):
                sub = sub[4:]
            if os.path.isdir(f'{self.output}/sub-{sub}/'):
                logging.warning(f'sub-{sub} preprocessing already completed!')
                if not overwrite:
                    logging.info(f'Skipping sub-{sub}')
                    continue
                else:
                    logging.warning(f'Re-preprocessing sub-{sub}, and overwriting results!')
            logging.info(f'Submitting Job {counter} of {len(self.subs)}')
            subprocess.run(f'bsub < {self.batch_dir}/sub-{sub}.sh', shell=True)
            counter += 1
            time.sleep(self.sleep_time)