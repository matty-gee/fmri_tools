function glm_design = glm_make_design(model, timing, behavior, glm_dir, verbose)

%------------------------------------------------------------------------------
% Make the design for a first-level GLM; also specifies some basic t-contrasts
% 
% 
% Arguments
% ---------
% model : cell array or str
%   specifies the model to build
% timing : table
%   contains onset information
% behavior : table
%   contains behavioral information: e.g., reaction-time, values for
%   parametric modulators etc
% glm_dir : str
%   where to output
% verbose : bool
% 
% 
% Returns
% -------
% glm_design : struct
% 
%
%
% Example model cell array: 
% >> model = {'all', 'decision', {'pov3d_angle','cos'}};
% 
% You can also enter a string for several common models: 
% 'lsa' : least squares all trialwise modeling with onset duration
% 'angle' : egocentric angle analysis from Tavares et al., 2015
% 'distance' : egocentric angle analysis from Tavares et al., 2015
% 'character' : 
% 'dimension' : 
%
%
%
% [By Matthew Schafer, github: @matty-gee; 2020ish] 
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
% helper functions/structures
%------------------------------------------------------------------------------


cond_struct = @(n,o,d,t) struct('name', n, 'onset', o, 'duration', d, 'trials', t, 'orth', 1,...
                                'pmod', struct('name', {}, 'param', {}, 'poly', {}, 'normalized', {}));  
event_durations = struct('onset', zeros(63,1), 'epoch', ones(63,1) * 12, 'decision', behavior.reaction_time);
event_weights   = zeros(63,1);
normalize = struct('none', @(X) X, 'z', @(X) zscore(X),...
                   'cos', @(X) cos(X), 'sin', @(X) sin(X),...
                   'cosd', @(X) cosd(X), 'sind', @(X) sind(X),...
                   'demean', @(X) X - mean(X));


%------------------------------------------------------------------------------
% define model 
%------------------------------------------------------------------------------


% if just a standard model name, create the model array now:
if ischar(model) 
    switch model
        case 'lsa'
            model  = {};
            for n_trial = 1:63
                model = [model; {sprintf('trial%02d', n_trial), 'onset', []}];
            end
        case 'character' 
            model = {};
            for n_char = 1:5
                model = [model; {sprintf('character0%d', n_char), 'decision', []}];
            end
        case 'dimensions'
            model = {'affil', 'decision', []; 'power', 'decision', []};
        case 'cos_angle'
            model = {'all', 'decision', {'pov_3d_angle','cos'}};
        case 'sin_angle'
            model = {'all', 'decision', {'pov_3d_angle','sin'}};
        case 'distance'
            model = {'all', 'decision', {'pov_dist','z'}};
        case 'polar_cos'
            model = {'all', 'decision', {'pov_3d_angle','cos'; 'pov_dist','z'}};
        case 'polar_sin'
            model = {'all', 'decision', {'pov_3d_angle','sin'; 'pov_dist','z'}};
        case 'polar_full'
            model = {'all', 'decision', {'pov_3d_angle','cos'; 'pov_3d_angle','sin'; 'pov_dist','z'}};
        case 'character_polar'
            model = {};
            for n_char = 1:5
                model = [model; {sprintf('character0%d', n_char), 'decision', {'pov_3d_angle', 'cos'; 'pov_dist', 'z'}}];
            end
    end
end

% if the model is character-wise w/ pmods, we need to specify each character condition
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
% define narrative condition first
%------------------------------------------------------------------------------


narrative = timing(strcmp(timing.trial_type, 'Narrative'), :);
regrs = {'narrative'};
glm_design.cond(1) = cond_struct('narrative', narrative.onset, narrative.duration, []);
n_conds = n_conds + 1; % add in narrative


% define the rest
decisions = timing(strcmp(timing.trial_type, 'Decision'), :);

for n_cond = 2:n_conds % excl. #1 b/c thats narrative condition

    [cond_trials, bold_event, cond_pmods] = model{n_cond-1, :}; % unpack the model


    %------------------------------------------------------------------------------
    % define event onsets
    %------------------------------------------------------------------------------


    if contains(cond_trials, 'trial') % single trial

        n_trial     = str2double(erase(cond_trials, 'trial'));
        cond_name   = sprintf('decision_%02d', n_trial);
        trials_incl = behavior.decision_num == n_trial;

    elseif strcmp(cond_trials, 'all')

        cond_name   = 'decisions';
        trials_incl = ones(63,1) == 1;

    elseif strcmp(cond_trials, 'affil') || strcmp(cond_trials, 'power')

        cond_name   = cond_trials;
        trials_incl = strcmp(behavior.dimension, cond_trials);

    elseif contains(cond_trials, 'character')

        n_char      = str2double(erase(cond_trials, 'character'));
        cond_name   = sprintf('character_%02d', n_char);
        trials_incl = behavior.char_role_num == n_char;

    end

    glm_design.cond(n_cond) = cond_struct(cond_name, decisions(trials_incl, :).onset, [], trials_incl);


    %------------------------------------------------------------------------------
    % assign durations [& optional parametric modulators] to events
    %------------------------------------------------------------------------------


    trials_incl = glm_design.cond(n_cond).trials;

    % duration
    glm_design.cond(n_cond).duration = event_durations.(bold_event)(trials_incl); 
    regrs = [regrs, glm_design.cond(n_cond).name];

    % pmods [optional]
    glm_design.cond(n_cond).orth = 0; % turn off orthogonalization - make sure uncorrelated...
    [n_cond_pmods, ~] = size(cond_pmods);
    for n_pmod = 1:n_cond_pmods

        [pmod_name, normlz] = cond_pmods{n_pmod, :};
        
        if ~strcmp(normlz, 'none')
           regr_name = [glm_design.cond(n_cond).name '*' normlz '(' pmod_name ')'];
        else
           regr_name = [glm_design.cond(n_cond).name '*' pmod_name]; 
        end
        regrs = [regrs, regr_name];
        param = behavior(trials_incl, pmod_name).Variables; 
        param = normalize.(normlz)(param);
        glm_design.cond(n_cond).pmod(n_pmod) = struct('name', pmod_name, 'param', param,...
                                                      'poly', 1, 'normalized', normlz);

    end
end

n_tps = length(vertcat(glm_design.cond.onset)); % number of timepoints


%------------------------------------------------------------------------------
% weight contrasts
%------------------------------------------------------------------------------


% weight individual regrs against baseline
n_regrs = length(regrs);
tcons = {};
for n_regr = 1:n_regrs

    tcon_weights = zeros(1, length(regrs));
    tcon_weights(n_regr) = 1;

    % positive 
    tcon_name = [regrs{n_regr} '+'];
    tcons = [tcons, tcon_name];
    glm_design.consess{length(tcons)}.tcon = struct('name', tcon_name, 'weights',  tcon_weights);

    % negative
    tcon_name = [regrs{n_regr} '-'];
    tcons = [tcons, tcon_name];
    glm_design.consess{length(tcons)}.tcon = struct('name', tcon_name, 'weights', -tcon_weights);
    
    % maybe add option to average related ones...?

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
    cond_ntrials = sum(horzcat(glm_design.cond(2:end).trials), 1);
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

f = filesep; % system-specific file separat

% save the model details
glm_design.model_array = model;
save([glm_dir f 'GLM.mat'], '-struct', 'glm_design')

% want to output a basic design matrix too, e.g., w canlab core tools

% text files w/ regrs & contrasts
writetable(cell2table(regrs), [glm_dir f 'GLM_regressors.csv'])
writetable(cell2table(tcons), [glm_dir f 'GLM_tcontrasts.csv'])

