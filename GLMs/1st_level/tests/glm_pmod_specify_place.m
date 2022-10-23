function specify_pmod_glm(sample, sub_id, pmod_name, nuisance, fwhm, hrf_name, hpf)

    % for place analysis, mainly

    addpath /hpc/packages/minerva-centos7/spm/spm12
    
    % need sub_id as a string
    if ~isstring(sub_id)
        if length(sub_id) == 1
            sub_id = sprintf('%02d', sub_id);
        else
            sub_id = num2str(sub_id);
        end
    end 
    
    % check that a viable model was entered
    possible_models = {'ego_angle', 'ego_distance', 'ego_angle_cf', 'ego_distance_cf'};
    if ~any(strcmp(possible_models, pmod_name))
        error('ERROR: specify parametric modulator: "ego_angle", "ego_distance", "ego_angle_cf" or "ego_distance_cf"')
    end
    
    % directories
    base_dir = '/sc/arion/projects/OlfMem/mgs/2D_place';    
    func_dir = [base_dir '/Samples/' sample '/spm_preprocessing/subs/' sub_id '/func'];
    glm_dir  = [base_dir '/Samples/' sample '/Pmods/' pmod_name '/' hrf_name... 
                '_' num2str(fwhm) '_hpf' num2str(hpf) '/subs/' sub_id];
           
    % check if model is already computed
    if isempty(spm_select('FPList', glm_dir, 'con_0001.nii')) 
        disp(['Preparing to run glm for ' sub_id])
    else
        disp(['Glm already computed for ' sub_id])
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % smoothed functional files
    if isempty(spm_select('FPList', func_dir, ['wau_func_smoothed' num2str(fwhm) '.nii']))    
        smooth_img(spm_select('ExtFPList', func_dir, 'wau_func.nii'), fwhm)
    end
    smoothed_imgs = cellstr(spm_select('ExtFPList', func_dir, ['wau_func_smoothed' num2str(fwhm) '.nii']));
    if isempty(smoothed_imgs), error("ERROR: Can't find func images"), end

    % nuisance regressors
    if strcmp(nuisance, 'rp')
        nuisance_txt = spm_select('FPList', func_dir, 'rp.txt');
        
    elseif strcmp(nuisance, 'rp_linear-trend')
        if isempty(spm_select('FPList', func_dir, 'rp_linear-trend.txt'))
            rp = spm_select('FPList', func_dir, 'rp.txt');
            cfds = [rp table((1:length(smoothed_imgs)).', 'VariableNames', {'linear_trend'})];
            writetable(cfds, [func_dir '/rp_linear-trend.txt'], 'WriteVariableNames', 0);
        end
        nuisance_txt = spm_select('FPList', func_dir, 'rp_linear-trend.txt');
    end
    if isempty(nuisance_txt), error("ERROR: Can't find nusiance text file"), end
    
    % get behavior for param modulator - also rt for pmod modeling 
    behav = readtable([base_dir '/Samples/' sample '/Behavior/Pmod_analyses/' sub_id '_pmods.xlsx']);
    behav = sortrows(behav, 'trial', 'ascend');
    
    % get onsets & durations - may vary a little between subjects
    timing = readtable([base_dir '/Samples/' sample '/Behavior/Timing/' sub_id '_timing.xlsx']);
    timing = sortrows(timing, 'onset', 'ascend');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% specify modeling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% microtime
    
    tr = [2; 1];
    mr = [36; 70];
    if strcmp(sample, 'Initial'), sample_num = 1;
    elseif strcmp(sample, 'Validation'), sample_num = 2;
    else, error("ERROR: specify a sample"),
    end
    pmod_model.tr = tr(sample_num);
    pmod_model.mr = mr(sample_num); % microtime resolution: if preprocessed imgs are slice time corrected, set to # of slices
    pmod_model.mo = mr(sample_num)/2; % microtime onset: temporal acq of ref slice (should prob. be middle slice)
    
    %% specify voxels to compute over
    
    pmod_model.mthresh = 0.5;
    pmod_model.mask_img = '';
   
    %% high-pass filter: filter out low frequency signal changes
    
    % -- default = 128s (1/128 Hz)
    % -- rule of thumb: threshlod at 2-3x average (or max) intervals (s) between predictor onsets (eg, https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;42334c51.1208)
    % -- Tavares 2015: 'below three cycles per time course" - 3/1500s --> 1/500s: this seems long..?
    
    pmod_model.hpf = hpf; % avg_interval = max(diff(pmod_model.decision_onsets))
    
    %% temporal autocorrelation regression: remove autocorrelated signal (ie prewhitening)
    % -- 'FAST' or 'AR(1)': FAST may be generally superior (Olszowy et al 2018)
    
    pmod_model.cvi = 'FAST';
    
    %% onset & offset modeling
    % -- if duration == 0, modeled as an event; otherwise as an epoch
    
    % decision onsets + rt
    pmod_model.decision_onsets = timing(strcmp(timing.trial_type, 'Decision'), :).onset;
    pmod_model.decision_durations = behav.rt; % 0 for onset, 12 for duration

    % narrative onsets + whole trial duration
    % -- excludes intro/outro image slides, intro & transition slides
    pmod_model.narrative_onsets = timing(strcmp(timing.trial_type, 'Narrative'), :).onset;
    pmod_model.narrative_durations = timing(strcmp(timing.trial_type, 'Narrative'), :).duration; 

    %% hemodynamic response function: model of neural events relationship to BOLD
    % -- [0 0]: canonical
    % -- [1 1]: temporal & dispersion derivatives
    % -- [1 0]: temporal derivative 
    
    % weighting beta images for contrasts
    % -- t-contrasts: should generally sum to 0 (unless contrasting against baseline)
    % -- modeling 'against baseline': baseline should prob. be left implicit to avoid over parameterization (Pernet 2014)
    % -- which/how many images to weight depends whether modeled partial derivatives of HRF (different number of images per condition)
    % ---- different arbitrary scaling between the images?
    % ---- if no reliable loading on canonical, difficult to interpret the partial derivatives
    % ---- might make most sense just to weight canonical, even if incl. derivatives
    
    if strcmp(hrf_name, 'canonHRF') == 1
        
        pmod_model.hrf = [0 0];
        
        pmod_model.consess{1}.tcon.name = 'decisions';
        pmod_model.consess{1}.tcon.weights = [1 0 0];
        pmod_model.consess{1}.tcon.sessrep = 'none';
        
        pmod_model.consess{2}.tcon.name = pmod_name;
        pmod_model.consess{2}.tcon.weights = [0 1 0];
        pmod_model.consess{2}.tcon.sessrep = 'none';
        
        pmod_model.consess{3}.tcon.name = 'narrative';
        pmod_model.consess{3}.tcon.weights = [0 0 1];
        pmod_model.consess{3}.tcon.sessrep = 'none';
        
    elseif strcmp(hrf_name, 'tempHRF') == 1
        
        % temporal derivative: allows peak response to vary by +/- 1s
        % partial derivative of canonical HRF wrt onset latency
        % -- 2 images produced per condition
        % -- [1 1] is average effect across those parameter estimates;
        % ---- most straightforward way to look at activation is canonical image
        % ---- (https://jiscmail.ac.uk/cgi-bin/wa.exe?A2=spm;af4afb49.1103)
        % for now: weight each separately so can be flexible on 2nd level
       
        pmod_model.hrf = [1 0]; 
       
        pmod_model.consess{1}.tcon.name = 'decisions_canon';
        pmod_model.consess{1}.tcon.weights = [1 0 0 0 0 0];
        pmod_model.consess{1}.tcon.sessrep = 'none';
        pmod_model.consess{2}.tcon.name = 'decisions_temp';
        pmod_model.consess{2}.tcon.weights = [0 1 0 0 0 0];
        pmod_model.consess{2}.tcon.sessrep = 'none';
        pmod_model.consess{3}.tcon.name = 'decisions';
        pmod_model.consess{3}.tcon.weights = [1 1 0 0 0 0];
        pmod_model.consess{3}.tcon.sessrep = 'none';
        
        pmod_model.consess{4}.tcon.name = [pmod_name '_canon'];
        pmod_model.consess{4}.tcon.weights = [0 0 1 0 0 0];
        pmod_model.consess{4}.tcon.sessrep = 'none';
        pmod_model.consess{5}.tcon.name = [pmod_name '_temp'];
        pmod_model.consess{5}.tcon.weights = [0 0 0 1 0 0];
        pmod_model.consess{5}.tcon.sessrep = 'none';
        pmod_model.consess{6}.tcon.name = pmod_name;
        pmod_model.consess{6}.tcon.weights = [0 0 1 1 0 0];
        pmod_model.consess{6}.tcon.sessrep = 'none';
        
        pmod_model.consess{7}.tcon.name = 'narrative_canon';
        pmod_model.consess{7}.tcon.weights = [0 0 0 0 1 0];
        pmod_model.consess{7}.tcon.sessrep = 'none';
        pmod_model.consess{8}.tcon.name = 'narrative_temp';
        pmod_model.consess{8}.tcon.weights = [0 0 0 0 0 1];  
        pmod_model.consess{8}.tcon.sessrep = 'none';
        pmod_model.consess{9}.tcon.name = 'narrative';
        pmod_model.consess{9}.tcon.weights = [0 0 0 0 1 1];  
        pmod_model.consess{9}.tcon.sessrep = 'none';
        
    elseif strcmp(hrf_name, 'tempdispHRF') == 1
        
        % partial derivatives of time and dispersion 
        
        pmod_model.hrf = [1 1];
        
        pmod_model.consess{1}.tcon.name = 'decisions_canon';
        pmod_model.consess{1}.tcon.weights = [1 0 0 0 0 0 0 0 0];
        pmod_model.consess{1}.tcon.sessrep = 'none';
        pmod_model.consess{2}.tcon.name = 'decisions_temp';
        pmod_model.consess{2}.tcon.weights = [0 1 0 0 0 0 0 0 0];
        pmod_model.consess{2}.tcon.sessrep = 'none';
        pmod_model.consess{3}.tcon.name = 'decisions_disp';
        pmod_model.consess{3}.tcon.weights = [0 0 1 0 0 0 0 0 0];
        pmod_model.consess{3}.tcon.sessrep = 'none';
        pmod_model.consess{4}.tcon.name = 'decisions';
        pmod_model.consess{4}.tcon.weights = [1 1 1 0 0 0 0 0 0];
        pmod_model.consess{4}.tcon.sessrep = 'none';
        
        pmod_model.consess{5}.tcon.name = [pmod_name '_canon'];
        pmod_model.consess{5}.tcon.weights = [0 0 0 1 0 0 0 0 0];
        pmod_model.consess{5}.tcon.sessrep = 'none';
        pmod_model.consess{6}.tcon.name = [pmod_name '_temp'];
        pmod_model.consess{6}.tcon.weights = [0 0 0 0 1 0 0 0 0];
        pmod_model.consess{6}.tcon.sessrep = 'none';
        pmod_model.consess{7}.tcon.name = [pmod_name '_disp'];
        pmod_model.consess{7}.tcon.weights = [0 0 0 0 0 1 0 0 0];
        pmod_model.consess{7}.tcon.sessrep = 'none';
        pmod_model.consess{8}.tcon.name = pmod_name;
        pmod_model.consess{8}.tcon.weights = [0 0 0 1 1 1 0 0 0];
        pmod_model.consess{8}.tcon.sessrep = 'none';
        
        pmod_model.consess{9}.tcon.name = 'narrative_canon';
        pmod_model.consess{9}.tcon.weights = [0 0 0 0 0 0 1 0 0];
        pmod_model.consess{9}.tcon.sessrep = 'none';
        pmod_model.consess{10}.tcon.name = 'narrative_temp';
        pmod_model.consess{10}.tcon.weights = [0 0 0 0 0 0 0 1 0];
        pmod_model.consess{10}.tcon.sessrep = 'none';
        pmod_model.consess{11}.tcon.name = 'narrative_disp';
        pmod_model.consess{11}.tcon.weights = [0 0 0 0 0 0 0 0 1];         
        pmod_model.consess{11}.tcon.sessrep = 'none';
        pmod_model.consess{12}.tcon.name = 'narrative';
        pmod_model.consess{12}.tcon.weights = [0 0 0 0 0 0 1 1 1];         
        pmod_model.consess{12}.tcon.sessrep = 'none';
        
    else
        
        error('ERROR: Something is off with the HRF/contrast weighting')
        
    end
      
    %% parametric modulator
    
    pmod_model.name = pmod_name;
    if strcmp(pmod_name, 'ego_angle') == 1
        pmod_model.param = cos(behav.ego_angles); % cosine similarity (angles in radians)
    elseif strcmp(pmod_name, 'ego_angle_cf') == 1
        pmod_model.param = cos(behav.ego_angles_cf); % counterfactual
    elseif strcmp(pmod_name, 'ego_distance') == 1
        pmod_model.param = zscore(behav.ego_distances); % z score euclidean distance
    elseif strcmp(pmod_name, 'ego_distance_cf') == 1
        pmod_model.param = zscore(behav.ego_distances_cf); % counterfactual
    end
    
    % polynomial expansion of the modulator:
    % -- if > 1, models non-linear relationships (2 = quadratic, 3 = cubic, etc)
    % -- im not sure how I would weight the resulting betas in a contrast...
    pmod_model.poly = 1;

    % parametric modulator orthogonalization:
    % -- if dont want to use, turn off spm orth & manualy demean 
    % -- spm serially orthogonalizes when multiple pmods - can be a weird inference
    pmod_model.orth = 1; 
%     if pmod_model.orth == 0
%        pmod_model.param = pmod_model.param - mean(pmod_model.param);
%     end
    
    %% add an explicit baseline?
    
    pmod_model.model_post_decisions = 0; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute glm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    compute_pmod_glm(smoothed_imgs, pmod_model, nuisance_txt, glm_dir)

end
    
