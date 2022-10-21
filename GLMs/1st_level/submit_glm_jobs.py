import json, os, glob, shutil, sys, time
import subprocess, logging

overwrite  = False
base_dir   = '/sc/arion/projects/k23'
batch_dir  = base_dir + '/code/GLMs/batch_dir'
preprc_dir = base_dir + '/derivatives_fieldmaps/fmriprep'
sub_dirs   = [sub for sub in glob.glob(preprc_dir + '/sub*') if 'html' not in sub] # exclude html files
sub_dirs.sort()

write_residuals = 0 
for glm in ['lsa', 'angle', 'distance', 'polar', 'affilpower_pmod']:
    for event in ['decision']:
        for nuisance in ['rp','rp+fd']:
            glm_dir = f'{base_dir}/GLMs_fieldmaps_{nuisance}/{glm}_{event}/subs'

            for sub_dir in sub_dirs:

                # create the job    
                sub_id = sub_dir.split('/')[-1]
                batch_script = f'{batch_dir}/{sub_id}_{glm}.sh'
                with open(batch_script, 'w') as f:
                    cookies = [ f'#!/bin/bash\n\n',
                                f'#BSUB -J {glm}_{sub_id}\n',
                                f'#BSUB -P acc_k23\n',
                                f'#BSUB -q premium\n',
                                f'#BSUB -n 2\n',
                                f'#BSUB -W 00:45\n',
                                f'#BSUB -R rusage[mem=8000]\n',
                                f'#BSUB -o {batch_dir}/batch_output/nodejob-{glm}-{sub_id}.out\n',
                                f'#BSUB -L /bin/bash\n\n',
                                f'ml matlab/R2021b\n\n', 
                                f'ml spm\n\n',
                                f'cd {base_dir}/code/GLMs/1st_level\n\n']
                    f.writelines(cookies)  
                    f.write(f"matlab -nodisplay -r \"glm_specify('{sub_id}', '{glm}', '{event}', '{nuisance}', {write_residuals});\"")

                # submit the job
                if os.path.isfile(f'{glm_dir}/{sub_id}/beta_0001.nii'):
                    logging.warning(f'{sub_id} already completed!')
                    if not overwrite:
                        logging.info(f'Skipping {sub_id}')
                        continue
                    else:
                        logging.warning(f'Re-runningg {sub_id}, and overwriting results!')
                logging.info(f'Submitting Job: {sub_id}')
                subprocess.run(f'bsub < {batch_script}', shell=True)
                time.sleep(600)
        # time.sleep(00)