import pandas as pd
import numpy as np
from nilearn.image import load_img
from sklearn.metrics import pairwise_distances
from scipy.stats import zscore

import json
from json import JSONEncoder
def load_json(filename):
    with open(filename, 'r') as f:
        return json.load(f)
        
def symm_mat_to_ut_vec(mat):
    vec_ut = mat[np.triu_indices(len(mat), k=1)]
    return vec_ut

# task details
task_details = pd.read_excel('task_details.xls') 
task_details.sort_values(by=['slide_num']) # make sure to sort
decision_details = task_details[task_details['trial_type']=='Decision']

# make control rdvs for rsa
time1_rdv = zscore(symm_mat_to_ut_vec(pairwise_distances(np.array(decision_details['onset']).reshape(-1,1), metric='euclidean')))
time2_rdv = zscore(symm_mat_to_ut_vec(pairwise_distances(np.array(decision_details['onset']).reshape(-1,1), metric='euclidean')**2))
time3_rdv = zscore(symm_mat_to_ut_vec(pairwise_distances(np.array(decision_details['onset']).reshape(-1,1), metric='euclidean')**3))
slide_rdv = zscore(symm_mat_to_ut_vec(pairwise_distances(np.array(decision_details['slide_num']).reshape(-1,1), metric='euclidean')))
scene_rdv = zscore(symm_mat_to_ut_vec(pairwise_distances(np.array(decision_details['scene_num']).reshape(-1,1), metric='euclidean')))
familiarity_rdv = zscore(symm_mat_to_ut_vec(pairwise_distances(np.array(decision_details['char_decision_num']).reshape(-1,1), metric='euclidean')))

def get_char_rdv(c):
    char_rdm = np.ones((63, 63))
    trial_nums = decision_details[decision_details['role_num'] == c]['decision_num'].values.astype(int) - 1
    for r in trial_nums:
        for c in trial_nums:
            char_rdm[c, r] = 0 
    return symm_mat_to_ut_vec(char_rdm)
char1_rdv = get_char_rdv(1)
char2_rdv = get_char_rdv(2)
char3_rdv = get_char_rdv(3)
char4_rdv = get_char_rdv(4)
char5_rdv = get_char_rdv(5)

def get_pred_df(beh_rdvs):
    
    decision_rdv = symm_mat_to_ut_vec(pairwise_distances(np.array(beh_rdvs['decisions']).reshape(-1,1))) 
    decision_rdv[decision_rdv > 1] = 1 # turn into catgorical 
    bp_rdv = np.array(beh_rdvs['bp_dv'])
    bp_rdv[bp_rdv > 1] = 1
        
    pred_df = pd.DataFrame({'place_2d': zscore(beh_rdvs['place_2d_dv']),\
                            'place_1d': zscore(beh_rdvs['place_1d_dv']),\
                            'place_affil': zscore(beh_rdvs['place_affil_dv']),\
                            'place_power': zscore(beh_rdvs['place_power_dv']),\
                            'place_positive': zscore(beh_rdvs['place_positive_dv']),\
                            'dist_ori': zscore(beh_rdvs['dist_ori_dv']),\
                            'dist_pov': zscore(beh_rdvs['dist_pov_dv']),\
                            'time1':time1_rdv,'time2':time2_rdv,'time3':time3_rdv,\
                            'slide': slide_rdv,'scene': scene_rdv,\
                            'familiarity': familiarity_rdv,\
                            'rt': zscore(beh_rdvs['rt_dv']),
                            'bp': bp_rdv,
                            'decisions': decision_rdv,\
                            'char1': char1_rdv,'char2': char2_rdv,'char3': char3_rdv,'char4': char4_rdv,'char5': char5_rdv,}) 

    return pred_df

# image loading & manipulation

def load_sub_imgs(sub_dir, img_name='beta_4d.nii', mask_name='mask.nii'):
    '''
        load data from img & brain mask
    '''
    img  = load_img(sub_dir + '/' + img_name)
    data  = img.get_fdata()
    affine_mat = img.affine
    vox_size = img.header.get_zooms() 
    try:
        brain_mask = load_img(sub_dir + '/' + mask_name).get_fdata()
    except:
        print('Missing a brain mask')
        brain_mask = None

    return data, brain_mask, affine_mat, vox_size

import nibabel as nib
# worse comes to worst just save as a 4d matrix np file and then convert to nifti later?
def save_as_nifti(brain_data, output_name, affine_mat, vox_size):
    
    brain_data = brain_data.astype('double')  # Convert the output into a precision format that can be used by other applications
    brain_data[np.isnan(brain_data)] = 0  # Exchange nans with zero to ensure compatibility with other applications
    brain_nii = nib.Nifti1Image(brain_data, affine_mat)  # create the volume image
    hdr = brain_nii.header  # get a handle of the .nii file's header
    if brain_data.ndim == 4: hdr.set_zooms((vox_size[0], vox_size[1], vox_size[2], 0))
    else: hdr.set_zooms((vox_size[0], vox_size[1], vox_size[2]))
    nib.save(brain_nii, output_name)  # Save the volume 