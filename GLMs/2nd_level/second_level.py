import pandas as pd
import numpy as np
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


#-------------------------------------------------------------------------
# 2nd level tests  
#-------------------------------------------------------------------------


def compute_ttest(fnames, design_matrix=None, fwhm=None, plot_design=False):
    '''
    '''
    print(f'Running t-test (n={len(fnames)}): fwhm={fwhm}')

    # design matrix
    if design_matrix is None: # assumes 1 sample t-test
        design_matrix = pd.DataFrame([1] * len(fnames), columns=['intercept']) 
    if plot_design:
        plot_design_matrix(design_matrix, rescale=False)
        plt.show()

    # estimate glm & make contrasts
    glm = SecondLevelModel(mask_img=None, smoothing_fwhm=fwhm, n_jobs=-1)
    glm = glm.fit(fnames, confounds=None, design_matrix=design_matrix)    
    
    # glm.generate_report  
    # # second_level_contrast=None,  output_type='z_score'
    # glm.compute_contrast(second_level_contrast=second_level_contrast, output_type=output_type)
    return glm

def threshold_second_level_img(con_img, 
                               mask_img=None, 
                               alpha=0.05, 
                               cluster_extent=0, 
                               two_sided=False, 
                               height_control='fdr'):
    
    ''' expects an image in z-scale '''

    print(f'Thresholding second level image: alpha={alpha}, cluster_extent={cluster_extent}, height_control={height_control}')
    con_img_thr, thr = threshold_stats_img(stat_img=con_img, 
                                         mask_img=mask_img,
                                         two_sided=two_sided, 
                                         alpha=alpha, 
                                         height_control=height_control,
                                         cluster_threshold=cluster_extent)
    return con_img_thr, thr

def compute_permutation_ttest(fnames, mask_img=None, fwhm=None, 
                              design_matrix=None, second_level_contrast=None, confound_df=None,
                              model_intercept=True, two_sided=False, threshold=None, 
                              n_perm=10000, tfce=False):

    ''' 
        Compute a permutation t-test on a list of first level contrast images

        Arguments
        ---------
        fnames: list of str
            list of paths to first level contrast images
        confound_df: pd.DataFrame
            dataframe with confounds, needs “subject_label” column
        mask_img: str
            path to mask image
        fwhm: float
            smoothing kernel in mm
        two_sided: bool
            whether to use two-sided test
        model_intercept: bool
            whether to model intercept
        n_perm: int
            number of permutations
        threshold: float
            p-scale cluster forming threshold
        tfce: bool
            whether to use threshold-free cluster enhancement

        Returns
        -------
        if threshold is None: negative logarithm of the voxel-level FWER-corrected p-values
        if threshold is not None: dictionary with keys (see: https://nilearn.github.io/dev/modules/generated/nilearn.glm.second_level.non_parametric_inference.html#nilearn.glm.second_level.non_parametric_inference)
    '''
    
    # https://nilearn.github.io/dev/modules/generated/nilearn.glm.second_level.non_parametric_inference.html

    print(f'Running nonparametric 1-sample t-test, n={len(fnames)}')

    if design_matrix is None: # assumes 1 sample t-test
        design_matrix = pd.DataFrame([1] * len(fnames), columns=['intercept']) 
    
    # produces a negative log pvalue image or an output dictionary if tfce=False and threshold is not None
    return non_parametric_inference(fnames, 
                                    design_matrix=design_matrix, 
                                    confounds=confound_df, # if used, needs “subject_label” column
                                    model_intercept=model_intercept, 
                                    first_level_contrast=None, 
                                    second_level_contrast=second_level_contrast, 
                                    mask=mask_img, 
                                    smoothing_fwhm=fwhm, 
                                    n_perm=n_perm, # number of 0s determines precision of the p-value
                                    two_sided_test=two_sided, 
                                    threshold=threshold, # p-scale cluster forming threshold
                                    tfce=tfce, # as described in orig. paper
                                    random_state=2022, n_jobs=-1, verbose=1) 

def threshold_neglog_img(neglog_img, alpha=0.05, cluster_threshold=0, vox_size=2.1):
    '''
        Threshold an image with negative logarithm (base 10) pvalues (larger is more significant)
    '''
    neglog_thresh = -np.log10(alpha)
    thr_img = nilearn.image.threshold_img(neglog_img, neglog_thresh, cluster_threshold=cluster_threshold)
    table = get_clusters_table(thr_img, stat_threshold=neglog_thresh, cluster_threshold=cluster_threshold) 
    if vox_size is not None:
        for i in range(len(table['Cluster Size (mm3)'])):
            c = table['Cluster Size (mm3)'][i]
            if isinstance(c, int): table['Cluster Size (mm3)'][i] = c / vox_size**3
        table = table.rename(columns={'Cluster Size (mm3)': 'Cluster Size (voxels)'})
    return thr_img, table


#-------------------------------------------------------------------------
# helpers
#-------------------------------------------------------------------------


def get_nifti_info(nii):
    if isinstance(nii, str):    
        nii = nib.load(nii)
    dims     = nii.get_fdata().shape
    vox_size = nii.header.get_zooms()[:3] # just get xyz
    affine   = nii.affine
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
