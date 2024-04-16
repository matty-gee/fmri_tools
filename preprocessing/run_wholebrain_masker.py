import os, glob, sys
base_dir = '/sc/arion/projects/OlfMem/mgs/2D_place'
sys.path.insert(0, f'{base_dir}/Code/toolbox') 
from images import save_as_nifti, get_nifti_info
import numpy as np
import nilearn

#---------------------------------------------------------------------------
# compute brain masks
#---------------------------------------------------------------------------

overwrite = True
base_dir = '/sc/arion/projects/OlfMem/mgs/2D_place'
for sample in ['Initial', 'Validation']:
    preprc_dir = os.path.join(base_dir, 'Samples', sample, 'spm_preprocessing', 'subs')
    subs = [d for d in glob.glob(f'{preprc_dir}/*') if not d.startswith('.')]
    for sub_dir in subs:
        sub_id = sub_dir.split('/')[-1]
        mask_fname = f'{sub_dir}/func/brain_mask.nii'
        if (os.path.exists(mask_fname)) and (not overwrite): 
            print(f'already computed brain mask for: {sub_id}', end='\r')
        else:
            print(f'masking {sub_id}', end='\r')

            brain_mask = nilearn.masking.compute_brain_mask(f'{sub_dir}/func/wau_func.nii', 
                                                            threshold=0.1, connected=True, opening=2, 
                                                            mask_type='whole-brain')
            dims, vox_size, affine = get_nifti_info(f'{sub_dir}/func/wau_func.nii')
            save_as_nifti(brain_mask.get_fdata(), mask_fname, affine, vox_size)