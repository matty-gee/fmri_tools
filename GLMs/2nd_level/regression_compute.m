function regression_compute(scans, cov, out_dir)

    % addpath /hpc/packages/minerva-centos7/spm/spm12
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');

    if ~exist(out_dir, 'dir'), mkdir(out_dir), end

    %% design matrix

    % scans
    batch{1}.spm.stats.factorial_design.dir = {out_dir};
    batch{1}.spm.stats.factorial_design.des.mreg.scans = scans;
    batch{1}.spm.stats.factorial_design.des.mreg.incint = 1; % include intercept
    batch{1}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
    
    % there are two cov related fields: factorial_design.des.mreg.mcov & factorial_design.cov
    % -- running in same way, with iCC = 1, yields same results...
    % -- but no iCFI option for interactions w/ factors

    % covariates in struct with fields:
    % -- cname: name of covariate
    % -- c: values
    % -- iCFI: covariate x factor interactions (creates an addtl regressor)
    % -- iCC: centering: makes intercept more interpretable
    batch{1}.spm.stats.factorial_design.cov = cov;

    % multiple covariates at once from a .mat file
    %  batch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});

    % masking
    batch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    batch{1}.spm.stats.factorial_design.masking.im = 1;
    batch{1}.spm.stats.factorial_design.masking.em = {''};

    batch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    batch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    batch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    %% estimation

    batch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    batch{2}.spm.stats.fmri_est.write_residuals = 0;
    batch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% contrasts
    
    batch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    % mean/intercept: effect controlling for covariates
    batch{3}.spm.stats.con.consess{1}.tcon.name    = 'intercept';
    batch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    
    % positive & negative effects
    for c = 1 : length(cov)
        
        con_weights = zeros(1, 1 + length(cov));
        
        % positive 
        con_weights(c+1) = 1;
        batch{3}.spm.stats.con.consess{2*c}.tcon.name    = [cov(c).cname '+'];
        batch{3}.spm.stats.con.consess{2*c}.tcon.weights = con_weights;

        % negative 
        con_weights(c+1) = -1;
        batch{3}.spm.stats.con.consess{2*c+1}.tcon.name    = [cov(c).cname '-'];
        batch{3}.spm.stats.con.consess{2*c+1}.tcon.weights = con_weights;
        
    end
    
    batch{3}.spm.stats.con.delete = 0;

    %% run

    spm_jobman('run', batch);