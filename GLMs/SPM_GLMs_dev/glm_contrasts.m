function glm_contrasts()

    batch{3}.spm.stats.con.spmmat                 = {[glm_dir f 'SPM.mat']};
    batch{3}.spm.stats.con.consess                = glm_model.consess;
    
    %% compute
    
    try
        spm_jobman('run', batch);
    catch
        disp('Problem running spm_jobman. Running in interactive mode.')
        spm_jobman('interactive', batch);
    end
    
    
        % if output residuals, concatenate residuals images
    % -- could restrict this to decision period, to make images smaller...
    if glm_model.write_residuals == 1
        res_3d = spm_select('FPList', glm_dir, '^Res_.*nii$');
        if isempty(res_3d), error('The residual images are missing'), end
        spm_file_merge(res_3d, 'Res_4d.nii');
        for r = 1 : length(res_3d), delete(res_3d(r, :)), end
    end

    % concatenate lsa images (lsa has 64 conditions: narrative (1) + 63 trials (2:64))
    if length(glm_model.cond) == 64
        beta_3d = spm_select('FPList', glm_dir, '^beta.*nii$');
        spm_file_merge(beta_3d(2 : 64, :), 'beta_4d.nii'); % skip first beta image: narrative condition
        tval_3d = spm_select('FPList', glm_dir, '^spmT.*nii$');
        spm_file_merge(tval_3d(1 : 63, :), 'spmT_4d.nii'); % dont have to skip first - only contrasts that were made
        for b = 1 : 63
            delete(beta_3d(b + 1, :))
            delete(tval_3d(b, :))
        end
    end