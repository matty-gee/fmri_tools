import pandas as pd
import matplotlib.pyplot as plt
import nibabel as nib

import nilearn
from nilearn.plotting import plot_design_matrix
from nilearn.image import get_data, math_img, threshold_img, binarize_img, new_img_like, smooth_img
from nilearn.input_data import NiftiMasker
from nilearn.masking import intersect_masks, compute_multi_brain_mask, compute_brain_mask, new_img_like
from nilearn.glm import cluster_level_inference, threshold_stats_img
from nilearn.glm.second_level import SecondLevelModel, non_parametric_inference
from nilearn.reporting import get_clusters_table

def compute_permutation_ttest(fnames, confound_df=None, mask_img=None, n_perm=1000, 
                              fwhm=None, cluster_threshold=None, two_sided=False, model_intercept=True, tfce=False):
    
    # https://nilearn.github.io/dev/modules/generated/nilearn.glm.second_level.non_parametric_inference.html
    
    print(f'{len(fnames)} subjects in a nonparametric 1-sample t-test')

    design_matrix = pd.DataFrame([1] * len(fnames), columns=['intercept']) 
    
    # produces a negative log pvalue image or an output dictionary if tfce=False and threshold is not None
    output = non_parametric_inference(fnames, 
                                      design_matrix=design_matrix, 
                                      confounds=confound_df, # if used, needs “subject_label” column
                                      model_intercept=model_intercept, 
                                      first_level_contrast=None, 
                                      second_level_contrast=None, 
                                      mask=mask_img, 
                                      smoothing_fwhm=fwhm, 
                                      n_perm=n_perm, # number of 0s determines precision of the p-value
                                      two_sided_test=two_sided, 
                                      threshold=cluster_threshold, # cluster forming threshold in p-scale
                                      tfce=tfce, # as described in the paper
                                      random_state=2022, n_jobs=-1, verbose=1) 

    return output

def compute_ttest(fnames, design_matrix=None, second_level_contrast=None,
                  fwhm=None, plot_design=True, output_type='z_score'):
    '''
    '''
    print(f'{len(fnames)} subjects in a t-test')

    # design matrix
    if design_matrix is None: # assume 1 sample t-test
        design_matrix = pd.DataFrame([1] * len(fnames), columns=['intercept']) 
    if plot_design:
        plot_design_matrix(design_matrix, rescale=False)
        plt.show()

    # estimate glm & make contrasts
    glm     = SecondLevelModel(mask_img=None, smoothing_fwhm=fwhm, n_jobs=-1)
    glm     = glm.fit(fnames, confounds=None, design_matrix=design_matrix)    
    con_img = glm.compute_contrast(second_level_contrast=second_level_contrast, output_type=output_type)
    
    # glm.generate_report
    
    return con_img

def threshold_second_level_img(con_img, mask_img=None, alpha=0.05, cluster_extent=0, two_sided=False, height_control='fdr'):
    '''
    '''
    con_img_thr, _ = threshold_stats_img(stat_img=con_img, 
                                         mask_img=mask_img,
                                         two_sided=two_sided, alpha=alpha, 
                                         height_control=height_control,
                                         cluster_threshold=cluster_extent)
    return con_img_thr

def get_nifti_info(nii):
    if isinstance(nii, str):    
        nii = nib.load(nii)
    dims = nii.get_fdata().shape
    vox_size = nii.header.get_zooms()[0:3] # just get xyz
    affine = nii.affine
    return dims, vox_size, affine 

def get_voxels_from_mask(mask_img, sub_img, resample_to_sub=False, standardize=False):
    '''
        mask_img: 3d nii (ideally already resampled to correct dims)
        sub_img: 4d nii
        returns: array of shape (time_points, voxels)
    '''
    if resample_to_sub:
        sub_dims, _, sub_affine = get_nifti_info(sub_img)
        masker = NiftiMasker(mask_img=mask_img, 
                             target_affine=sub_affine, target_shape=sub_dims[0:3],
                             standardize=standardize)
    else:
        masker = NiftiMasker(mask_img=mask_img, standardize=standardize)
    return masker.fit_transform(sub_img)

def get_perm_pval(coef, perm_coefs, test='2-sided'):
    '''2-sided, greater or lesser
    '''
    p_greater = np.sum(perm_coefs >= coef) / len(perm_coefs)
    p_less  = np.sum(perm_coefs <= coef) / len(perm_coefs)
    p_2sided  = min(p_greater, p_less)
    if test == '2-sided':
        pval = p_2sided
    elif test == 'greater':
        pval = p_greater
    elif test == 'less':
        pval = p_less
    return pval
