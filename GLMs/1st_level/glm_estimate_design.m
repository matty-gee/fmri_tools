function glm_estimate_design(spm_mat, write_residuals)

    spm('defaults', 'FMRI')

    batch{1}.spm.stats.fmri_est.spmmat(1)         = spm_mat;
    batch{1}.spm.stats.fmri_est.write_residuals   = write_residuals;
    batch{1}.spm.stats.fmri_est.method.Classical  = 1;

    try
        spm_jobman('run', batch);
    catch
        disp('Problem running spm_jobman. Running in interactive mode.')
        spm_jobman('interactive', batch);
    end