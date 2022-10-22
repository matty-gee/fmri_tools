function glm_specify_design(func_dir, model, write_residuals, output, verbose)

%------------------------------------------------------------------------------
% 
% 
% Arguments
% ---------
% func_dir
%
% model
%
% write_residuals
% output
% verbose
%
%
% Accepted model names:
%
%
%
% Example model array:
% use ';' to separate diff conditions and diff pmods within condition
% eg: a model with one condition but 2 pmods: 
% model = {cond_trials, bold_event [optional],...
%         {cond_pmod1 [optional], cond_pmod1_normalization [optional];...
%          pmod2 [optional], pmod2_normalization [optional];}}
%
% eg: a model with two conditions and 1 pmod each: 
% model = {cond1_trials, bold_event [optional],...
%         {cond1_pmod1 [optional], cond1_pmod1_normalization [optional]};... 
%         {cond2_trials, bold_event [optional],...
%         {cond2_pmod1 [optional], cond2_pmod1_normalization [optional]};... 
%
%
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% unpack the model
%------------------------------------------------------------------------------   


smooth = 1; % most models will require smoothing

% if just a standard model name, create the model array now:
if ischar(model) 
    switch model
        case 'lsa'
            smooth = 0;
            model  = {};
            for n_trial = 1:63
                model = [model; {sprintf('trial%02d', n_trial), 'onset', []}];
            end
        case 'character' 
            model = {};
            for n_char = 1:5
                model = [model; {sprintf('character0%d', n_char), 'decision',[]}];
            end            
        case 'angle'
            model = {'all', 'decision', {'pov3d_angle','cos'}};
        case 'distance'
            model = {'all', 'decision', {'pov_distance','z'}};
        case 'dimensions'
            model = {'affil', 'decision', []; 'power', 'decision', []};
    end
end

% if character-wise model w/ pmods, specify each character condition
if strcmp(model{1}, 'character') 
    tmp = {};
    for n_char = 1:5
        tmp = [tmp; {sprintf('character0%d', n_char), model{2}, model{3}}];
    end
    model = tmp;
end

[n_conds, ~] = size(model);
[n_pmods, ~] = size(vertcat(model{:, 3}));

%------------------------------------------------------------------------------
% directories
%------------------------------------------------------------------------------

addpath /sc/arion/projects/k23/code/fmri_utilities % add circ_mean here
addpath /hpc/packages/minerva-centos7/spm/spm12

f = filesep; % system-specific 
 
% base_dir = [f 'sc' f 'arion' f 'projects' f 'k23'];
% func_dir = [base_dir f 'derivatives_fieldmaps' f 'fmriprep' f sub_id f 'func']; 
% beh_dir  = [base_dir f 'behavior']; 
% glm_dir  = [base_dir f 'GLMs_fieldmaps_' cfds f glm_name f bold_event f 'subs' f sub_id]; % for the output
% if ~exist(glm_dir, 'dir'), mkdir(glm_dir), end
% disp(['Preparing to run GLM for ' sub_id '...'])

%------------------------------------------------------------------------------
% project/analysis - can i read this stuff from the nifti instead? tr & mr?
%------------------------------------------------------------------------------


glm.tr  = 1; 
glm.mr  = 70; % microtime resolution: number of slices collected
glm.mo  = glm.mr/2; 
glm.hrf = [0 0]; 
glm.hpf = 128; % default = 128s (1/128 Hz) [or: max(diff(pmod_model.decision_onsets))]
glm.cvi = 'FAST'; % prewhitening to remove autocorrelated signal 


%------------------------------------------------------------------------------
% images, masking & nuisance regrs
%------------------------------------------------------------------------------


% funcional images
func_name = 'wau_func'; 
func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_name '.nii$']));
% img_unzip(func_dir, func_fname)

% if smooth
%     img_smooth(spm_select('ExtFPList', func_dir, [func_name '.nii']), 6)
%     func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_fname '_smoothed6.nii']));  
% end

% masking
glm.mthresh  = 0.50; % default=0.8; -Inf allows for explicit masks
glm.mask_img = '';
% img_unzip(func_dir, '^s.*brain_mask')
% glm_model.mthresh  = -Inf; % default=0.8; -Inf allows for explicit masks
% glm_model.mask_img = spm_select('FPList', func_dir, '^s.*brain_mask.nii$');

% nuisance regrs
nuisance_txt = [func_dir f 'rp.txt']; 
% nuisance_txt = cfd_fmriprep_nuisance_txt(func_dir, cfds);

%------------------------------------------------------------------------------
% events: behavior, onsets & durations
%------------------------------------------------------------------------------


%behav = readtable([beh_dir f 'Pmod_analyses' f erase(sub_id, 'sub-P') '_pmods.xlsx']);
behav = '/Users/matty_gee/Dropbox/Projects/social_navigation_task/data/example_subject/snt_18001_behavior.xlsx';
behav = readtable(behav);
behav = sortrows(behav, 'decision_num', 'ascend');

% timing = readtable([beh_dir f 'Timing' f erase(sub_id, 'sub-P') '_timing.xlsx']); % change to hardcode if single timing file
timing = '/Users/matty_gee/Dropbox/Projects/social_navigation_task/data/example_subject/snt_18001_timing.xlsx';
timing = readtable(timing);
timing = sortrows(timing, 'onset', 'ascend');


%------------------------------------------------------------------------------s
% helper functions/structures
%------------------------------------------------------------------------------


cond_struct = @(n,o,d,t) struct('name', n, 'onset', o, 'duration', d, 'trials', t, 'orth', 1,...
                                'pmod', struct('name', {}, 'param', {}, 'poly', {}, 'normalized', {}));  
event_durations = struct('onset', zeros(63,1), 'epoch', ones(63,1) * 12, 'decision', behav.reaction_time);
event_weights   = zeros(63,1);
normalize = struct('none', @(X) X, 'z', @(X) zscore(X), 'cos', @(X) cos(X), 'cosd', @(X) cosd(X));


%------------------------------------------------------------------------------
% define narrative condition first
%------------------------------------------------------------------------------


narrative = timing(strcmp(timing.trial_type, 'Narrative'), :);
regrs = {'narrative'};
glm.cond(1) = cond_struct('narrative', narrative.onset, narrative.duration, []);
n_conds = n_conds + 1; % add in narrative


% define the rest
for n_cond = 2:n_conds % x isnt informative in anyway

    [cond_trials, bold_event, cond_pmods] = model{n_cond-1, :}; % unpack the model


    %------------------------------------------------------------------------------
    % define event onsets
    %------------------------------------------------------------------------------


    if contains(cond_trials, 'trial') % single trial

        n_trial     = str2double(erase(cond_trials, 'trial'));
        cond_name   = sprintf('decision_%02d', n_trial);
        trials_incl = behav.decision_num == n_trial;

    elseif strcmp(cond_trials, 'all')

        cond_name   = 'decisions';
        trials_incl = ones(63,1) == 1;

    elseif strcmp(cond_trials, 'affil') || strcmp(cond_trials, 'power')

        cond_name   = cond_trials;
        trials_incl = strcmp(behav.dimension, cond_trials);

    elseif contains(cond_trials, 'character')

        n_char      = str2double(erase(cond_trials, 'character'));
        cond_name   = sprintf('character_%02d', n_char);
        trials_incl = behav.char_role_num == n_char;

    end

    glm.cond(n_cond) = cond_struct(cond_name, behav(trials_incl, :).onset, [], trials_incl);


    %------------------------------------------------------------------------------
    % assign durations [& optional parametric modulators] to events
    %------------------------------------------------------------------------------

    trials_incl = glm.cond(n_cond).trials;

    % duration
    glm.cond(n_cond).duration = event_durations.(bold_event)(trials_incl); 

    % pmods
    regrs = [regrs, glm.cond(n_cond).name]; 
    [n_cond_pmods, ~] = size(cond_pmods);
    for n_pmod = 1:n_cond_pmods % optional

        [pmod_name, normlz] = cond_pmods{n_pmod, :};
        regrs = [regrs, [glm.cond(n_cond).name '*' pmod_name]];
        param = behav(trials_incl, pmod_name).Variables; 

        % multiple pmods: turn off orthogonalization, demean
        if length(cond_pmods) > 1
            glm.cond(n_cond).orth = 0;
            if contains(pmod_name, 'angle'), param = param - circ_mean(param); 
            else,                            param = param - mean(param); 
            end
        end

        param = normalize.(normlz)(param); % normalize: after demeaning b/c cosine isnt true "normalization"
        glm.cond(n_cond).pmod(n_pmod) = struct('name', pmod_name, 'param', param,...
                                               'poly', 1, 'normalized', normlz);

    end
end
n_tps = length(vertcat(glm.cond.onset));


%------------------------------------------------------------------------------
% weight contrasts
%------------------------------------------------------------------------------


% weight individual regrs (1) against baseline
n_regrs = length(regrs);
tcons = {};
for n_regr = 1:n_regrs

    tcon_weights = zeros(1, length(regrs));
    tcon_weights(n_regr) = 1;

    % positive 
    tcon_name = [regrs{n_regr} '+'];
    tcons = [tcons, tcon_name];
    glm.consess{length(tcons)}.tcon = struct('name', tcon_name, 'weights',  tcon_weights);

    % negative
    tcon_name = [regrs{n_regr} '-'];
    tcons = [tcons, tcon_name];
    glm.consess{length(tcons)}.tcon = struct('name', tcon_name, 'weights', -tcon_weights);

end
n_tcons = length(tcons);


%------------------------------------------------------------------------------
% do some checks each time
%------------------------------------------------------------------------------


check_sum = @(to_sum, exp_sum) sum(to_sum) == exp_sum;

% check number of events in diff conditions (excls narrative)
if n_conds-1 ~= 63 % dont do lsa
    if n_conds-1 == 1
        cond_ntrials_exp = 63;
    elseif n_conds-1 == 2 % assuming a power/affil breakdown
        cond_ntrials_exp = [30 30];
    elseif n_conds-1 == 5 % assuming characterwise 
        cond_ntrials_exp = [12 12 12 12 12];
    end
    cond_ntrials = sum(horzcat(glm.cond(2:end).trials), 1);
    if check_sum(cond_ntrials == cond_ntrials_exp, n_conds-1) == 0
        error('Number of trials per condition is off')
    end
end

% check number of regressors (incls narrative)
if check_sum((n_pmods + n_conds), n_regrs) == 0
    error('Number of regressors is off')
end

if verbose 
    fprintf('%d timepoints\n', n_tps)
    fprintf('%d conditions with %d total param. mod.\n', n_conds, n_pmods)
    fprintf('%d regressors: %s\n', n_regrs, sprintf('%s; ', regrs{:}))
    fprintf('%d contrasts: %s\n', n_tcons, sprintf('%s; ', tcons{:}))
end


%------------------------------------------------------------------------------
% output stuff
%------------------------------------------------------------------------------


if output

    glm.model_array = model;
    save([glm_dir f 'GLM_' model_name '.mat'], '-struct', 'glm')

    % want to output a basic design matrix too, e.g., w canlab core tools

    % text files w/ regrs & contrasts
    writetable(cell2table(regrs), [glm_dir f 'GLM_' model_name '_regrs.csv'])
    writetable(cell2table(tcons), [glm_dir f 'GLM_' model_name '_regrs.csv'])
end


%------------------------------------------------------------------------------
% compute glm
%------------------------------------------------------------------------------


% clean up glm object so can pass glm.cond directly into spm
% glm = rmfield(glm, 'model_array'); % wont get used by glm_compute.m
glm.cond = rmfield(glm.cond, 'trials');
for n_cond = 1:length(glm.cond)
    glm.cond(n_cond).pmod = rmfield(glm.cond(n_cond).pmod, 'normalized');
end

glm.model_post_decisions = 0;
glm.write_residuals      = write_residuals; % write residuals (1) or not (0) [eg, 1 if want to do func conn after]
% glm_compute(func_imgs, glm, nuisance_txt, glm_dir)


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% dev
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

% TODO: for more flexible contrast weighting - 
% f = @(r, p) (contains(r,'*') == p) * 1; % are regrs conds (p=0) or pmods (p=1)
% pmods_inc = cumsum(f(regrs,1)) .* f(regrs,1); % number pmods
% conds_inc = cumsum(f(regrs,0)) .* f(regrs,0);

end