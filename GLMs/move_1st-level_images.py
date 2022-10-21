import glob
import nibabel as nib
import os
import shutil

glm_dict = {'angle_decision':['con_0003.nii'], 
            'distance_decision':['con_0003.nii'],
            'polar_decision': ['con_0003.nii', 'con_0004.nii'], # pmods: angle, decision
            'lsa_decision': ['beta_4d.nii', 'mask.nii']} 

for glm, img_fnames in glm_dict.items():
    glm_dir = '/sc/arion/projects/k23/GLMs_fieldmaps_rp+fd+csf/' + glm
    if not os.path.exists(glm_dir + '/images'): os.makedirs(glm_dir + '/images')
    for sub_dir in glob.glob(glm_dir + '/subs/*'):
        print(sub_dir)
        sub_id = sub_dir.split('/')[-1]
        for img_fname in img_fnames:
            try:    shutil.copy(sub_dir + '/' + img_fname, 
                                glm_dir + '/images/' + sub_id + '_' + img_fname) 
            except: print(f'error with {sub_dir}/{img_fname}')