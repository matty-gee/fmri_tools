function glm_specify(sub_id, glm_name, bold_event, cfds, write_residuals)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specifies first-level GLM
    %
    % Arguments
    % ---------
    % sub_id : str
    % glm_name : str 
    %   Available GLMs defined under 'supported GLMs'
    % bold_event : str
    %   Event to model: 'onset', 'decision', or 'epoch'
    % cfds : str
    %   Nuisance variables to include, separated by '+' (e.g., 'rp+fd')
    %   See: cfd_fmriprep_nuisance_txt.m
    % write_residuals: boolean
    %   Set = 1 if want residuals 
    %
    % EXAMPLE: subject 18001, angle parametric modulation, over the decision period (ie, rt), with realignment
    % parameters and not outputting residuals
    % >> glm_specify('sub-P18001', 'angle', 'decision', 'rp', 0)
    % 
    % [By Matthew Schafer, github: @matty-gee; 2020ish] 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% supported GLMs:    
    % - 'lsa': each trial (63) is own condition 
    % - 'simple': options & narrative conditions
    % - 'affilpower_cond': affil & power conditions
    % - 'affilpower_cond_pmod': affil & power conditions w/ pmod
    % - 'affilpower_pmod': decision condition w/ affil & power pmods
    % - 'angle', 'distance', 'polar': decision condition w/ pmod(s)

    % directories

    addpath /sc/arion/projects/k23/code/utilities
    addpath /hpc/packages/minerva-centos7/spm/spm12

    f = filesep; % system-specific 
    base_dir = [f 'sc' f 'arion' f 'projects' f 'k23'];
    func_dir = [base_dir f 'derivatives_fieldmaps' f 'fmriprep' f sub_id f 'func']; 
    beh_dir  = [base_dir f 'behavior']; 
    glm_dir  = [base_dir f 'GLMs_fieldmaps_' cfds f glm_name f bold_event f 'subs' f sub_id]; % for the output
    if ~exist(glm_dir, 'dir'), mkdir(glm_dir), end
    
    disp(['Preparing to run GLM for ' sub_id '...'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% project/analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % scan timing, hrf 
    glm_model.tr  = 1; 
    glm_model.mr  = 70; % microtime resolution: number of slices collected
    glm_model.mo  = glm_model.mr/2; 
    glm_model.hrf = [0 0]; 
    
    % filtering  
    glm_model.hpf = 128; % default = 128s (1/128 Hz) [or: max(diff(pmod_model.decision_onsets))]
    glm_model.cvi = 'FAST'; % prewhitening to remove autocorrelated signal 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% subject data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % functional images (unzip & smooth if needed)
    func_fname = '^s.*preproc_bold'; 
    img_unzip(func_dir, func_fname)
    
    if strcmp(glm_name, 'lsa') 
        func_imgs = cellstr(spm_select('ExtFPList', func_dir, '^s.*preproc_bold.nii$')); % use unsmoothed images
        glm_model.hrf = [0 0]; % use canonical
    else
        img_smooth(spm_select('ExtFPList', func_dir, [func_fname '.nii$']), 6)
        func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_fname '_smoothed6.nii$']));  
    end

    % masking
    img_unzip(func_dir, '^s.*brain_mask')  
    glm_model.mthresh  = -Inf; % default=0.8; -Inf allows for explicit masks
    glm_model.mask_img = spm_select('FPList', func_dir, '^s.*brain_mask.nii$');
    
    % if need to intersect multiple masks... use spm's imcalc(m1 * m2 * m3 * m4)
    % 

    % nuisance regressors
    nuisance_txt = cfd_fmriprep_nuisance_txt(func_dir, cfds);
    %     if contains(cfds, 'acc') % dont use hpf option if using aCompCor since cosine regressors will be included in nuisance file
    %         glm_model.hpf = % but no way to turn off spm's hpf?
    %     end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% events: behavior, onsets & durations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    behav  = readtable([beh_dir f 'Pmod_analyses' f erase(sub_id, 'sub-P') '_pmods.xlsx']);
    behav  = sortrows(behav, 'trial', 'ascend');
    timing = readtable([beh_dir f 'Timing' f erase(sub_id, 'sub-P') '_timing.xlsx']); % change to hardcode if single timing file
    timing = sortrows(timing, 'onset', 'ascend');

    % timing: onsets & durations
    narrative_timing = timing(strcmp(timing.trial_type, 'Narrative'), :);
    decision_onsets  = timing(strcmp(timing.trial_type, 'Decision'), :).onset;
    affil_onsets     = timing(strcmp(timing.dimension, 'affil'), :).onset;
    power_onsets     = timing(strcmp(timing.dimension, 'power'), :).onset;
    
    switch bold_event
        
        case 'onset' % onset of decision trial
            durations       = zeros(63,1);
            affil_durations = zeros(30,1);
            power_durations = zeros(30,1);
        case 'epoch' % epoch of decision trial
            durations       = ones(63,1)*12;
            affil_durations = ones(30,1)*12;
            power_durations = ones(30,1)*12;
        case 'decision' % reaction time of decision
            durations       = behav.rt;
            affil_durations = behav(strcmp(behav.dimension, 'affil'), :).rt;
            power_durations = behav(strcmp(behav.dimension, 'power'), :).rt;
            
    end
    
    % if want to filter out missed trials, do something like this: 
    % missed_trials_mask = durations > 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% define events & contrasts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    glm_model.cond(1) = spm_cond_struct('narrative', narrative_timing.onset, narrative_timing.duration); % 1st condition: narrative
    
    switch glm_name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trialwise conditions (e.g., for MVPA)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 'lsa' % least-squares all [e.g., nilearn tutorials on lsa v. lss in beta-series modeling]
            
            for t = 1 : 63
                glm_model.cond(1+t) = spm_cond_struct(['trial_' num2str(t, '%02d')], decision_onsets(t), durations(t)); 
                glm_model.consess{t}.tcon = struct('name', ['trial_' num2str(t, '%02d')], 'weights', (0:64 == t) * 1); % [narr, trial_1, ..., trial_63]
            end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % decision condition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case {'simple', 'angle', 'angle_cf', 'distance', 'distance_cf', 'polar', 'affilpower_pmod'} 
        
            glm_model.cond(2) = spm_cond_struct('decision', decision_onsets, durations);
            glm_model.consess{1}.tcon = struct('name', 'narrative > decision', 'weights', [1 -1]);
            glm_model.consess{2}.tcon = struct('name', 'decision > narrative', 'weights', [-1 1]);

            switch glm_name 
                
                % 1 pmod
                case {'angle','angle_cf','distance','distance_cf'}

                    switch glm_name
                        case 'angle'
                            pmod = cos(behav.ego_angles); % cosine similarity; angles in radians -> use cos function
                        case 'angle_cf' % cf = counterfactual
                            pmod = cos(behav.ego_angles_cf); 
                        case 'distance'
                            pmod = zscore(behav.ego_distances); % z-scored
                        case 'distance_cf'
                            pmod = zscore(behav.ego_distances_cf); 
                    end
                    
                    glm_model.cond(2).pmod    = struct('name', glm_name, 'param', pmod, 'poly', 1);
                    glm_model.consess{3}.tcon = struct('name', [glm_name '+'], 'weights', [0 0 1]); % [narr, dec, dec*pmod]
                    glm_model.consess{4}.tcon = struct('name', [glm_name '-'], 'weights', [0 0 -1]);
            
                % 2 pmods
                case {'polar', 'affilpower_pmod'}

                    if strcmp(glm_name, 'polar') % angle & distance
                        names = {'angle', 'distance'};
                        pmod1 = cos(behav.ego_angles);
                        pmod2 = zscore(behav.ego_distances);
                    elseif strcmp(glm_name, 'affilpower_pmod') % affil & power
                        names = {'affiliation', 'power'};
                        pmod1 = zscore(behav.affil_coord); 
                        pmod2 = zscore(behav.power_coord);  
                    end

                    glm_model.cond(2).orth = 0; % no orthogonalization 
                    glm_model.cond(2).pmod(1) = struct('name', names{1}, 'param', pmod1-mean(pmod1), 'poly', 1); % demean
                    glm_model.cond(2).pmod(2) = struct('name', names{2}, 'param', pmod2-mean(pmod2), 'poly', 1);

                    glm_model.consess{3}.tcon = struct('name', [names{1} '+'],            'weights', [0 0 1]); % [narr, dec, dec*pmod1, dec*pmod2]
                    glm_model.consess{4}.tcon = struct('name', [names{1} '-'],            'weights', [0 0 -1]);
                    glm_model.consess{5}.tcon = struct('name', [names{2} '+'],            'weights', [0 0 0 1]);
                    glm_model.consess{6}.tcon = struct('name', [names{2} '-'],            'weights', [0 0 0 -1]); 
                    glm_model.consess{7}.tcon = struct('name', [names{1} ' > ' names{2}], 'weights', [0 0 1 -1]); 
                    glm_model.consess{8}.tcon = struct('name', [names{2} ' > ' names{1}], 'weights', [0 0 -1 1]); 

            end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % affil & power decision conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case {'affilpower_cond', 'affilpower_cond_pmod'}
            
            glm_model.cond(2) = spm_cond_struct('affiliation', affil_onsets, affil_durations);
            glm_model.cond(3) = spm_cond_struct('power',       power_onsets, power_durations);
            
            if strcmp(glm_name, 'affilpower_cond')
                weights = [0 1 2];     % [narr, affil, power]
            elseif strcmp(glm_name, 'affilpower_cond_pmod')
                weights = [0 0 1 0 2]; % [narr, affil, affil*pmod, power, power*pmod]
                glm_model.cond(2).pmod = struct('name', glm_name, 'param', zscore(behav.affil_coord), 'poly', 1);
                glm_model.cond(3).pmod = struct('name', glm_name, 'param', zscore(behav.power_coord), 'poly', 1);  
            end

            glm_model.consess{1}.tcon = struct('name', 'affiliation',         'weights', (weights == 1) * 1);                        % [0 1 0],  [0 0 1]
            glm_model.consess{2}.tcon = struct('name', 'power',               'weights', (weights == 2) * 1);                        % [0 0 1],  [0 0 0 0 1]
            glm_model.consess{3}.tcon = struct('name', 'affiliation > power', 'weights', (weights == 1) * 1  + (weights == 2) * -1); % [0 1 -1], [0 0 1 0 -1]
            glm_model.consess{4}.tcon = struct('name', 'power > affiliation', 'weights', (weights == 1) * -1 + (weights == 2) * 1);  % [0 -1 1], [0 0 -1 0 1]

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % character conditions 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 'character_angle'
       
            n_con = 1;
            c_weights = [0 1 0 2 0 3 0 4 0 5 0];
            p_weights = [0 0 1 0 2 0 3 0 4 0 5];
            
            for c = 1 : 5

                character = behav(behav.role == c, :);
                glm_model.cond(1+c) = spm_cond_struct(['character_' num2str(c, '%02d')], character.onset, character.rt);
                glm_model.cond(1+c).pmod = struct('name', 'angle', 'param', character.ego_angles, 'poly', 1);

                glm_model.consess{n_con}.tcon   = struct('name', ['character_' num2str(c, '%02d')],           'weights', (c_weights == c) * 1);
                glm_model.consess{n_con+1}.tcon = struct('name', ['character_' num2str(c, '%02d') '_angle+'], 'weights', (p_weights == c) * 1);

                n_con = n_con + 2;
                
            end

            % average
            glm_model.consess{11}.tcon = struct('name', 'characters',        'weights', (c_weights ~= 0) * 1);
            glm_model.consess{11}.tcon = struct('name', 'characters_angle+', 'weights', (p_weights ~= 0) * 1);

    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute glm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    glm_model.model_post_decisions = 0;
    glm_model.write_residuals      = write_residuals; % write residuals (1) or not (0) [eg, 1 if want to do func conn after]
    glm_compute(func_imgs, glm_model, nuisance_txt, glm_dir)
