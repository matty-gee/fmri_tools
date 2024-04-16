function ttest_ind_compute(scans1, scans2, out_dir)

    addpath /hpc/packages/minerva-centos7/spm/spm12
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');
    
    if ~exist(out_dir, 'dir'), mkdir(out_dir), end

    % design matrix
    batch{1}.spm.stats.factorial_design.dir = {out_dir};
    batch{1}.spm.stats.factorial_design.des.t2.scans1 = scans1;
    batch{1}.spm.stats.factorial_design.des.t2.scans2 = scans2;
    batch{1}.spm.stats.factorial_design.des.t2.dept = 0;
    batch{1}.spm.stats.factorial_design.des.t2.variance = 1;
    batch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
    batch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
    batch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    batch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    batch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    batch{1}.spm.stats.factorial_design.masking.im = 1;
    batch{1}.spm.stats.factorial_design.masking.em = {''};
    batch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    batch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    batch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % estimation
    batch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    batch{2}.spm.stats.fmri_est.write_residuals = 0;
    batch{2}.spm.stats.fmri_est.method.Classical = 1;

    % contrasts
    batch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    batch{3}.spm.stats.con.consess{1}.fcon.name = 'different';
    batch{3}.spm.stats.con.consess{1}.fcon.weights = [-1 1];
    batch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    
    batch{3}.spm.stats.con.consess{2}.tcon.name = 'greater';
    batch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1];
    batch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    batch{3}.spm.stats.con.consess{3}.tcon.name = 'lesser';
    batch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1];
    batch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';    

    batch{3}.spm.stats.con.delete = 0;

%     % report w/ FDR correction
%     batch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     batch{4}.spm.stats.results.conspec.titlestr = '';
%     batch{4}.spm.stats.results.conspec.contrasts = 1;
%     batch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
%         % if want voxel-wise FDR, edit the spm_defaults file to read defaults.stats.topoFDR = 0;
%     batch{4}.spm.stats.results.conspec.thresh = 0.05;
%     batch{4}.spm.stats.results.conspec.extent = 0;
%     batch{4}.spm.stats.results.conspec.conjunction = 1;
%     batch{4}.spm.stats.results.conspec.mask.none = 1;
%     batch{4}.spm.stats.results.units = 1;
%     batch{4}.spm.stats.results.export{1}.ps = true;
%     
    % run
    spm_jobman('run', batch);
    
end