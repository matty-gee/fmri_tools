import json, os, glob, shutil, sys   
import subprocess, logging
import time

overwrite  = False
counter    = 1
base_dir   = '/sc/arion/projects/k23/VBM'
sub_dirs   = [sub for sub in glob.glob(base_dir + '/sub*')] 

for sub_dir in sub_dirs:

    sub_id = sub_dir.split('/')[-1]
    sub_batch_script = f'{base_dir}/batch_dir/{sub_id}.sh'
    anat_dir = sub_dir + '/anat'
    
    # create the job
    with open(sub_batch_script, 'w') as f:
        cookies = [ f'#!/bin/bash\n\n',
                    f'#BSUB -J {sub_id}\n',
                    f'#BSUB -P acc_guLab\n', 
                    f'#BSUB -q private\n', 
                    f'#BSUB -n 2\n',
                    f'#BSUB -W 00:30\n',
                    f'#BSUB -R rusage[mem=8000]\n',
                    f'#BSUB -o {base_dir}/batch_dir/batch_output/nodejob-{sub_id}.out\n',
                    f'#BSUB -L /bin/bash\n\n',
                    f'ml matlab/R2021b\n\n', 
                    f'ml spm\n\n',
                    f'cd {base_dir}/code/VBM\n\n']
        f.writelines(cookies)  
        command = f"matlab -nodisplay -r \"vbm_run('{anat_dir}');\""
        f.write(command)

    # submit the job
    if os.path.isfile(f'{anat_dir}/mwp1{sub_id}_T1w_smoothed8.nii'):
        logging.warning(f'{sub_id} already completed!')
        if not overwrite:
            logging.info(f'Skipping {sub_id}')
            continue
        else:
            logging.warning(f'Re-preprocessing {sub_id} for VBM, and overwriting results!')
    logging.info(f'Submitting Job {counter} of {len(sub_dirs)}: {sub_id}')
    subprocess.run(f'bsub < {sub_batch_script}', shell=True)
    counter += 1
    time.sleep(60)
time.sleep(300)
