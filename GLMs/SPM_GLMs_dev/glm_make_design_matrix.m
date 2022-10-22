function glm_make_design_matrix(func_imgs, glm_model, nuisance_txt, glm_dir)

    spm('defaults', 'FMRI')

    % io
    batch{1}.spm.stats.fmri_spec.dir              = {glm_dir}; 
    batch{1}.spm.stats.fmri_spec.sess.scans       = func_imgs; 

    % timing 
    batch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    batch{1}.spm.stats.fmri_spec.timing.RT        = glm_model.tr;
    batch{1}.spm.stats.fmri_spec.timing.fmri_t    = glm_model.mr; % microtime resolution: how many time-bins to use per volume to build regressors; if stc, set to number of slices
    batch{1}.spm.stats.fmri_spec.timing.fmri_t0   = glm_model.mo; % microtime onset: if stc, match ref slice (probably middle slice)

    % condtions (incl. parametric modulation)
    for c = 1 : length(glm_model.cond)        
        batch{1}.spm.stats.fmri_spec.sess.cond(c) = glm_model.cond(c);
    end

    % nuisance regressors
    batch{1}.spm.stats.fmri_spec.sess.multi_reg   = {nuisance_txt};

    % filtering, masking etc
    batch{1}.spm.stats.fmri_spec.sess.hpf         = glm_model.hpf; % default = 128s (1/128 Hz); rule of thumb: threshold at 2-3x average (or max) intervals (s) between predictor onsets
    batch{1}.spm.stats.fmri_spec.bases.hrf.derivs = glm_model.hrf; 
    batch{1}.spm.stats.fmri_spec.mthresh          = glm_model.mthresh; % '-Inf': all voxels
    batch{1}.spm.stats.fmri_spec.mask             = {glm_model.mask_img}; % explicit mask
    batch{1}.spm.stats.fmri_spec.cvi              = glm_model.cvi; % pre-whitening of serial correlations: default='AR(1)' or 'FAST'; FAST might be better at removing temporally autocorrelated BOLD signal (Olszowy et al 2018)

    spm_jobman('run', batch);