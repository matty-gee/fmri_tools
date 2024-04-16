import json, os, glob, shutil, sys, time
import subprocess, logging
import pandas as pd

base_dir   = '/sc/arion/projects/OlfMem/mgs/2D_place'

# for job files
batch_dir  = f'{base_dir}/Code/extract_timeseries/batch_dir'
for dir_ in ['jobs', 'out']:
    if not os.path.exists(f'{batch_dir}/{dir_}'):
        os.makedirs(f'{batch_dir}/{dir_}')

# diff options for job submission
projects  = ['acc_guLab', 'acc_k23']
queues    = ['private', 'premium']
submit_to = 1

# specify some glm stuff
write_residuals = 0 
overwrite       = False

# subjects
participants = pd.read_excel('participants_info.xlsx', usecols=['sub_id']).values.astype(int)

for sample in ['Validation']:

    # pick a subject
    preprc_dir = f'{base_dir}/Samples/{sample}/spm_preprocessing/subs/'
    sub_dirs   = [sub for sub in glob.glob(f'{preprc_dir}/*') if 'html' not in sub] # exclude html files
    sub_dirs.sort()
    
    for sub_dir in sub_dirs:
        sub_id = sub_dir.split('/')[-1]
        
        # create the job
        batch_sh = f'{batch_dir}/jobs/{sub_id}.sh'
        with open(batch_sh, 'w') as f:
            cookies = [ f'#!/bin/bash\n\n',
                        f'#BSUB -J {sub_id}\n',
                        f'#BSUB -P {projects[submit_to]}\n',
                        f'#BSUB -q {queues[submit_to]}\n',
                        f'#BSUB -n 2\n',
                        f'#BSUB -W 00:20\n',
                        f'#BSUB -R rusage[mem=8000]\n',
                        f'#BSUB -o {batch_dir}/out/{sub_id}.out\n',
                        f'#BSUB -L /bin/bash\n\n',
                        f'ml matlab/R2021b\n\n', 
                        f'ml spm\n\n',
                        f'cd {base_dir}/Code/GLMs/\n\n']
            f.writelines(cookies)
            f.write(f"matlab -nodisplay -r \"extract_timeseries('{sub_dir}');\"") 

        logging.info(f'Submitting Job: {sub_id}')
        subprocess.run(f'bsub < {batch_sh}', shell=True)
        time.sleep(20)
