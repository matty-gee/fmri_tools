function glm_design_matrix_png(glm_dir)
    
    %visualize
    load([glm_dir filesep 'SPM.mat'])
     
    % individal regressors
    for c = 1:length(SPM.Sess.U)
    	spm_DesRep('fMRIDesMtx', SPM, 1, 2) % should be decisions trial
        saveas(gcf, [glm_dir filesep 'Design_condition_' SPM.Sess.U(c).name{1} '.png'])
    end
    
    % overall design matrix
    spm_DesRep('DesMtx', SPM.xX)
    saveas(gcf, [glm_dir filesep 'Design_matrix.png'])