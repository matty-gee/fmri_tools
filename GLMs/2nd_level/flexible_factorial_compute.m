function flexible_factorial_compute(subjects, cov, out_dir)

    % specify design matrix flexibly
    % best for anovas where dont want to test all possible main & interaction fx
    % 3 stages: specify factors, add scans, build design matrix

    % addpath /hpc/packages/minerva-centos7/spm/spm12
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');

    if ~exist(out_dir, 'dir'), mkdir(out_dir), end
    batch{1}.spm.stats.factorial_design.dir = {out_dir};
    
    %% specify design matrix
    
    % specify factor(s): name, in/dependence, equal variance
    % -- for group, should prob. assume independence and unequal variance
    % -- main fx: produces # of regressors = factor levels 
    % -- interaction fx:  produces # of regressors = factor levels * condition
    batch{1}.spm.stats.factorial_design.des.fblock.fac.name     = 'group'; 
    batch{1}.spm.stats.factorial_design.des.fblock.fac.dept     = 0; % 0=independence
    batch{1}.spm.stats.factorial_design.des.fblock.fac.variance = 1; % 1=unequal variance
    
    % specify subjects: single measurement per
    for s = 1 : length(subjects)
        batch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {subjects(s).scans}; 
        batch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = subjects(s).conds; % specify which factor & level each image belongs to
    end
    
    % main effects & interactions
    batch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
%     batch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = '<UNDEFINED>'; % interactions between factors
    
    % specify covariates with a struct with fields:
    % -- cname: name
    % -- c: values
    % -- iCFI: covariate x factor interactions (creates addtl regressors)
    % -- iCC: centering: makes intercept more interpretable
    batch{1}.spm.stats.factorial_design.cov = cov;
    
    % masking, etc
    batch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1; % no threshold masking
%     batch{1}.spm.stats.factorial_design.masking.tm.tma.athresh = 100; % absolute threshold
%     batch{1}.spm.stats.factorial_design.masking.tm.tmr.rthresh = 0.8; % relative threshold
    batch{1}.spm.stats.factorial_design.masking.im             = 1; % implicit mask
    batch{1}.spm.stats.factorial_design.masking.em             = {''}; % explicit mask
    batch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
    batch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    batch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

    %% estimate glm
    
    batch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File',... 
                                                    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
                                                    substruct('.','spmmat'));
    batch{2}.spm.stats.fmri_est.write_residuals  = 0;
    batch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% weight regressors for contrasts
    
    batch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
                                                substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}),... 
                                                substruct('.','spmmat'));
    
    % group factor
    
    % -- groups by themselves
    batch{3}.spm.stats.con.consess{1}.tcon.name    = 'Gr1+';
    batch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0];
    batch{3}.spm.stats.con.consess{2}.tcon.name    = 'Gr1-';
    batch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 0];
    
    batch{3}.spm.stats.con.consess{3}.tcon.name    = 'Gr2+';
    batch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
    batch{3}.spm.stats.con.consess{4}.tcon.name    = 'Gr2-';
    batch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    
    % -- contrasted
    batch{3}.spm.stats.con.consess{5}.tcon.name    = 'Gr1 > Gr2';
    batch{3}.spm.stats.con.consess{5}.tcon.weights = [1 -1];
    batch{3}.spm.stats.con.consess{6}.tcon.name    = 'Gr2 > Gr1';
    batch{3}.spm.stats.con.consess{6}.tcon.weights = [-1 1];
    
    % -- grand means
    batch{3}.spm.stats.con.consess{7}.tcon.name    = 'Gr Avg+';
    batch{3}.spm.stats.con.consess{7}.tcon.weights = [1 1];
    batch{3}.spm.stats.con.consess{8}.tcon.name    = 'Gr Avg-';
    batch{3}.spm.stats.con.consess{8}.tcon.weights = [-1 -1];
    
    % positive & negative effects of covariates
    
    n_cons = length(batch{3}.spm.stats.con.consess); % contrast count
    reg    = 2; % which regressor
    for c = 1 : length(cov)
        
        % weight the parameter estimates to create linear contrasts
        
        if cov(c).iCFI == 1
            
            reg = reg + 1;
            
            n_cons           = n_cons + 1;
            con_weights      = zeros(1, reg); 
            con_weights(reg) = 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname '+'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;

            n_cons           = n_cons + 1;
            con_weights(reg) = -1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname '-'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;
            
        elseif cov(c).iCFI == 2
            
            % interaction w/ a factor creates 2 regressors
            
            reg1 = reg + 1; % regressor 1
            reg2 = reg + 2; % regressor 2

            % covariate x group 1
            con_weights       = zeros(1, reg1); 
            con_weights(reg1) = 1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr1+'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;

            con_weights(reg1) = -1; % flip sign
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr1-'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;
            
            % covariate x group 2
            con_weights       = zeros(1, reg2); 
            
            con_weights(reg2) = 1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr2+'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;

            con_weights(reg2) = -1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr2-'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;

            % contrast covariate x group fx 
            con_weights       = zeros(1, reg2); 
            con_weights(reg1) = 1;
            con_weights(reg2) = -1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr1 > ' cov(c).cname 'xGr2'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;

            con_weights(reg1) = -1;
            con_weights(reg2) = 1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr2 > ' cov(c).cname 'xGr1'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;
            
            % avg. covariate x group effects
            con_weights       = zeros(1, reg2); 
            con_weights(reg1) = 1;
            con_weights(reg2) = 1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr Avg+'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;

            con_weights(reg1) = -1;
            con_weights(reg2) = -1;
            n_cons            = n_cons + 1;
            batch{3}.spm.stats.con.consess{n_cons}.tcon.name    = [cov(c).cname 'xGr Avg-'];
            batch{3}.spm.stats.con.consess{n_cons}.tcon.weights = con_weights;
            
            reg = reg + 2; % add 2 regressors
        end

    end
    
    %% run

    spm_jobman('run', batch);
   