function glm_compute(func_imgs, glm_design, nuisance_txt, glm_dir)

%------------------------------------------------------------------------------------
% Computes first-level GLM
%
% Arguments
% ---------
% func_imgs : cell array
%   functional images to model
% glm_design : struct
%   matlab structure
% nuisance_txt : str
% glm_dir : str
%   
% [By Matthew Schafer, github: @matty-gee; 2020ish] 
%------------------------------------------------------------------------------------


% minerva paths, no harm no foul
addpath /sc/arion/projects/k23/code/matlab_utilities
addpath /hpc/packages/minerva-centos7/spm/spm12

f = filesep; % system specific


%------------------------------------------------------------------------------------
% error checking etc
%------------------------------------------------------------------------------------


% check the nifti length
% if (length(func_imgs) ~= 1570) && (length(func_imgs) ~= 784 && (length(func_imgs) ~= 848))
%     error(['ERROR: Wrong number of func. volumes: ' num2str(length(func_imgs)) ' instead of 784, 848 or 1570.'])
% end

run_it = 1;

% delete any old spm.mats if incompleted
if (isfile([glm_dir f 'SPM.mat'])) && (~isfile([glm_dir f 'con_0001.nii'])) % (~isfile([glm_dir f 'RPV.nii']))
    disp('Unfinished old estimation. Deleting old SPM.mat')
    delete([glm_dir f 'SPM.mat'])

% skip if already completed
elseif (isfile([glm_dir f 'SPM.mat'])) && (isfile([glm_dir f 'con_0001.nii']))
    disp('Estimation appears to have been completed already. Exiting.')
    run_it = 0;
end

spm('defaults', 'FMRI')
disp('Running GLM...')
if glm_design.write_residuals == 1, disp('Outputting residuals'), end

if run_it == 1
    
    %------------------------------------------------------------------------------------
    % define design matrix
    %------------------------------------------------------------------------------------


    % io
    batch{1}.spm.stats.fmri_spec.dir              = {glm_dir}; 
    batch{1}.spm.stats.fmri_spec.sess.scans       = func_imgs; 

    % timing 
    batch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    batch{1}.spm.stats.fmri_spec.timing.RT        = glm_design.tr;
    batch{1}.spm.stats.fmri_spec.timing.fmri_t    = glm_design.mr; % microtime resolution: how many time-bins to use per volume to build regressors; if stc, set to number of slices
    batch{1}.spm.stats.fmri_spec.timing.fmri_t0   = glm_design.mo; % microtime onset: if stc, match ref slice (probably middle slice)

    % condtions (incl. parametric modulation)
    for c = 1 : length(glm_design.cond)        
        batch{1}.spm.stats.fmri_spec.sess.cond(c) = glm_design.cond(c);
    end

    % nuisance regressors
    batch{1}.spm.stats.fmri_spec.sess.multi_reg   = {nuisance_txt};

    % filtering, masking etc
    batch{1}.spm.stats.fmri_spec.sess.hpf         = glm_design.hpf; % default = 128s (1/128 Hz); rule of thumb: threshold at 2-3x average (or max) intervals (s) between predictor onsets
    batch{1}.spm.stats.fmri_spec.bases.hrf.derivs = glm_design.hrf; 
    batch{1}.spm.stats.fmri_spec.mthresh          = glm_design.mthresh; % '-Inf': all voxels
    batch{1}.spm.stats.fmri_spec.mask             = {glm_design.mask_img}; % explicit mask
    batch{1}.spm.stats.fmri_spec.cvi              = glm_design.cvi; % pre-whitening of serial correlations: default='AR(1)' or 'FAST'; FAST might be better at removing temporally autocorrelated BOLD signal (Olszowy et al 2018)


    %------------------------------------------------------------------------------------
    % estimate glm
    %------------------------------------------------------------------------------------


    batch{2}.spm.stats.fmri_est.spmmat(1)         = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    batch{2}.spm.stats.fmri_est.write_residuals   = glm_design.write_residuals;
    batch{2}.spm.stats.fmri_est.method.Classical  = 1;


    %------------------------------------------------------------------------------------
    % weight contrasts
    %------------------------------------------------------------------------------------


    batch{3}.spm.stats.con.spmmat                 = {[glm_dir f 'SPM.mat']};
    batch{3}.spm.stats.con.consess                = glm_design.consess;


    %------------------------------------------------------------------------------------
    % compute
    %------------------------------------------------------------------------------------


    try
        spm_jobman('run', batch);
    catch
        disp('Problem running spm_jobman. Running in interactive mode.')
        spm_jobman('interactive', batch);
    end


    %------------------------------------------------------------------------------------
    % clean up files
    %------------------------------------------------------------------------------------

    % if outputting residuals, concatenate 3D residuals images into 4D image
    if glm_design.write_residuals == 1
        res_3d = spm_select('FPList', glm_dir, '^Res_.*nii$');
        if isempty(res_3d), error('The residual images are missing'), end
        spm_file_merge(res_3d, 'Res_4d.nii');
        for r = 1 : length(res_3d), delete(res_3d(r, :)), end
    end

    % if lsa, concatenate 3D lsa images into 4D image (lsa has 64 conditions: narrative (1) + 63 trials (2:64))
    if length(glm_design.cond) == 64

        beta_3d = spm_select('FPList', glm_dir, '^beta.*nii$');
        spm_file_merge(beta_3d(2 : 64, :), 'beta_4d.nii'); % skip first beta image: narrative condition

        tval_3d = spm_select('FPList', glm_dir, '^spmT.*nii$');
        spm_file_merge(tval_3d(2 : 64, :), 'spmT_4d.nii'); 

        % clean up
        for b = 2 : 64
            delete(beta_3d(b, :)) 
            delete(tval_3d(b, :))
        end
    end

end
% rename beta & con files appropiately... subid_betaname.nii, etc...