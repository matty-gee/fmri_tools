def run_sl_rsa(sub_dir, gm_mask, radius=3, shape='ball', min_prop=.10):

    import os
    import time
    import glob
    import warnings
    import pandas
    import nilearn
    import numpy as np    
    warnings.filterwarnings('ignore')

    from nilearn.image import load_img
    from sklearn.linear_model import HuberRegressor
    from sklearn.feature_selection import VarianceThreshold
    import brainiak.searchlight.searchlight
    from utils_sl_minerva import load_json, get_pred_df, load_sub_imgs, symm_mat_to_ut_vec, save_as_nifti

    def rsa_sl(brain_data, sl_mask, myrad, bcvar):
        '''
            inputs:
                brain_data: the voxels
                sl_mask: what voxels to run searchlight over
                
                need these variables for the brainiak searchlight object: 
                myrad: radius of searchlight
                bcvar: 
        '''    
        
        # neural data -> searchlight 
        data_4D = brain_data[0] 
        data_sl = data_4D.reshape(sl_mask.shape[0] * sl_mask.shape[1] * sl_mask.shape[2], brain_data[0].shape[3]).T
        
        try: 
            # get correlation distances
            data_sl = VarianceThreshold().fit_transform(data_sl) # variance threshold to remove all 0 or nan cols
            neural_rdv = 1 - symm_mat_to_ut_vec(np.arctanh(np.corrcoef(data_sl))) # fisher-z transform the correlations

            # fit model 
            regression = HuberRegressor(epsilon=1.75, alpha=0.0001, max_iter=50)
            regression.fit(predictor_rdvs, neural_rdv)
            
            # keep the coefficient for the predictor of interest
            if isinstance(pred_pos, slice):
                coef = np.mean(regression.coef_[pred_pos]) # get an average
            else: 
                coef = regression.coef_[pred_pos]
        except: 
            coef = np.nan
            
        return coef

    controls = ['bp', 'rt', 'char1', 'char2', 'char3', 'char4', 'char5', 
                'slide', 'scene', 'familiarity',  'time1', 'time2', 'time3']
    predictors = ['place_2d'] + controls
    pred_pos = 0 # slice(3,8)

    # searchlight properties
    sl_shapes = {'cube': brainiak.searchlight.searchlight.Cube,
                 'ball': brainiak.searchlight.searchlight.Ball,
                 'diamond': brainiak.searchlight.searchlight.Diamond}
    max_blk_edge = 10 # number of sls to run concurrently
    pool_size = 1 # number of cores (usually 1)sub_dir + '/' + sub_mask_name

    # define directories, file names etc 
    if 'Original' in sub_dir:
        sample = 'Original' 
        img_name = 'beta_4d_resampled.nii'
        sub_mask_name = 'mask_resampled.nii'
    else:
        sample = 'Replication' 
        img_name = 'beta_4d.nii'
        sub_mask_name = 'mask.nii'   
    sample_dir = '../../Samples/' + sample    

    # what kind of mask
    if gm_mask is None:
        sl_dir = f'{sample_dir}/GLMs_FAST_mask50/Searchlights/place_2d_{radius}-voxel-{shape}_{min_prop}'
    else:
        if isinstance(gm_mask, float):
            mask_name = 'thr' + str(gm_mask)           
        elif isinstance(gm_mask, str):
            if gm_mask == 'All-subjs_gm_thr25.nii':
                mask_name = 'gm-all-thr25'
            elif gm_mask == 'All-subjs_gm_thr50.nii':
                mask_name = 'gm-all-thr50'
        sl_dir = f'{sample_dir}/GLMs_FAST_mask50/Searchlights/place_2d_{radius}-voxel-{shape}_{min_prop}_{mask_name}'

    # get subject's data  
    sub_id = sub_dir.split('/')[-1].split('_')[0]
    sl_fname = os.path.join(sl_dir, ('%s_sl.nii' % (sub_id)))   
    beh_rdv_dict = load_json(sample_dir + '/Behavior/Behavioral_distances.json')        
    predictor_rdvs = get_pred_df(beh_rdv_dict[sub_id])[predictors]    
    brain_img, _, affine_mat, vox_size = load_sub_imgs(sub_dir, img_name=img_name, mask_name=sub_mask_name)

    # how to mask?
    if gm_mask is None: # just with the 1st level spm mask
        brain_mask = nilearn.image.load_img(sub_dir + '/' + sub_mask_name)
    else: # intersect the subject level mask with a GM mask 
        if isinstance(gm_mask, str): # use predefined mask    
            brain_mask = nilearn.image.binarize_img(sub_dir + '/' + sub_mask_name, threshold=0, mask_img=load_img('../../Masks/' + gm_mask)) 
        elif isinstance(gm_mask, float): # use predefined threshold to compute gm mask from sub mask
            # minerva doesnt recognize 'masking' rn: maybe an old version of nilearn?
            brain_mask = nilearn.masking.compute_brain_mask(sub_dir + '/' + sub_mask_name, mask_type='gm', threshold=gm_mask, connected=False)      
    brain_mask = brain_mask.get_fdata() # need data as matrix to mask searchlight

    # build & distribute sl to run
    print('sub_id', sub_id)
    print('func data shape', brain_img.shape) 
    print('mask shape', brain_mask.shape, '\n') 
    print(shape, 'radius', radius)
    begin_time = time.time()
    sl = brainiak.searchlight.searchlight.Searchlight(sl_rad=radius, shape=sl_shapes[shape], 
                                                      min_active_voxels_proportion=min_prop,
                                                      max_blk_edge=max_blk_edge) 
    sl.distribute([brain_img], brain_mask) # searchlight only in the brain mask                             
    sl_result = sl.run_searchlight(rsa_sl, pool_size=pool_size) 
    save_as_nifti(sl_result, sl_fname, affine_mat, vox_size)
    print('Total searchlight duration (including start up time): %.2f' % (time.time() - begin_time))
