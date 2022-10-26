import json, os, glob, shutil, sys, time
import subprocess, logging

base_dir   = '/sc/arion/projects/OlfMem/mgs/2D_place'

# for job files
batch_dir  = f'{base_dir}/Code/GLMs/batch_dir'
for dir_ in ['jobs', 'out']:
    if not os.path.exists(f'{batch_dir}/{dir_}'):
        os.makedirs(f'{batch_dir}/{dir_}')

# specify some glm stuff
write_residuals = 0 
overwrite       = False

for glm in ['angle', 'distance']:
    for sample in ['Initial', 'Validation']:

        # pick a subject
        glm_dir    = f'{base_dir}/Samples/{sample}/GLMs/{glm}'
        preprc_dir = f'{base_dir}/Samples/{sample}/spm_preprocessing/subs/'
        sub_dirs   = [sub for sub in glob.glob(f'{preprc_dir}/*') if 'html' not in sub] # exclude html files
        sub_dirs.sort()
        
        for sub_dir in sub_dirs:

            # create the job    
            sub_id   = sub_dir.split('/')[-1]
            batch_sh = f'{batch_dir}/jobs/{sub_id}_{glm}.sh'
            with open(batch_sh, 'w') as f:
                cookies = [ f'#!/bin/bash\n\n',
                            f'#BSUB -J {sub_id}_{glm}\n',
                            f'#BSUB -P acc_k23\n',
                            f'#BSUB -q premium\n',
                            f'#BSUB -n 2\n',
                            f'#BSUB -W 00:45\n',
                            f'#BSUB -R rusage[mem=8000]\n',
                            f'#BSUB -o {batch_dir}/out/{sub_id}_{glm}.out\n',
                            f'#BSUB -L /bin/bash\n\n',
                            f'ml matlab/R2021b\n\n', 
                            f'ml spm\n\n',
                            f'cd {base_dir}/Code/GLMs/\n\n']
                f.writelines(cookies)  
                # args: sub_dir, model, glm_name, verbose, debug
                f.write(f"matlab -nodisplay -r \"glm_specify_design_spmprep('{sub_dir}', '{glm}', '{glm}', {write_residuals}, 0);\"")

            # submit the job
            print(f'{glm_dir}/subs/{sub_id}')
            if os.path.isfile(f'{glm_dir}/subs/{sub_id}/beta_0001.nii'):
                logging.warning(f'{sub_id} already completed!')
                if not overwrite:
                    logging.info(f'Skipping {sub_id}')
                    continue
                else:
                    logging.warning(f'Re-runningg {sub_id}, and overwriting results!')
            logging.info(f'Submitting Job: {sub_id}')
            subprocess.run(f'bsub < {batch_sh}', shell=True)
            time.sleep(20)
