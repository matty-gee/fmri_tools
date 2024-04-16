# move images on minerva to cleanup a bit 

import os, glob, shutil
import nibabel as nib

base_dir = '/sc/arion/projects/OlfMem/mgs/2D_place/'
out_dir  = f'{base_dir}/lsa_images'
glm_dict = {'lsa_decision': ['beta_0001.nii', 'beta_4d.nii', 'spmT_0001.nii', 'spmT_4d.nii', 'mask.nii'],
            'lsa_onset': ['beta_0001.nii', 'beta_4d.nii', 'spmT_0001.nii', 'spmT_4d.nii', 'mask.nii']}

for sample in ['Initial', 'Validation']:
    for glm, img_fnames in glm_dict.items():
        glm_dir = f'{base_dir}/Samples/{sample}/GLMs/{glm}'
        for sub_dir in glob.glob(f'{glm_dir}/subs/*'):
            print(sub_dir)
            sub_id = sub_dir.split('/')[-1]
            for img_fname in img_fnames:
                try:    
                    shutil.copy(f'{sub_dir}/{img_fname}', f'{out_dir}/{glm}_128hpf_{sub_id}_{img_fname}') 
                except: 
                    print(f'error with {sub_dir}/{img_fname}')