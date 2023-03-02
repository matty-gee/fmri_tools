%------------------------------------------------------------------------------------------
% run diff models through glm_specify to test whether we get correct shapes
%------------------------------------------------------------------------------------------

clear 

sub_dir = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/Scans/spm_preprocessing/subs/22010';

% trials, bold_evemt, {pmod, pmod_normalization}
model1 = {'all', 'decision', []}; % 1 cond, no pmods
model2 = {'all', 'decision', {'pov3d_angle','cos'}}; % 1 cond, 1 pmod
model3 = {'all', 'decision', {'pov3d_angle','cos'; 'pov_distance','z'}}; % 1 cond, 2 pmods
model4 = {'affil', 'decision', {'affil_mean','z'};... % 2 conds w/ own pmods
          'power', 'decision', {'power_mean','z'}};  
model5 = {'character', 'decision', {'pov3d_angle','cos'}}; % 5 conds ('characters'), w/ same 1 pmod
model6 = {'affil', 'decision', {'affil_mean','z'; 'power_mean','z'};... % 2 conds w/ own pmods
          'power', 'decision', {'power_mean','z'}};

% weird models to test edges
% - should throw an error
model7 = {'all', 'decision', {'affil_mean','z'; 'power_mean','z'};... 
          'power', 'decision', {'power_mean','z'}};
      
% models = [{model1}; {model2}; {model3}; {model4}; {model5}; {model6}; {'lsa'}; {'character'};...
%           {'angle'}; {'distance'}; {'dimensions'};];

models = [{'polar_full'}; {'character_polar'}];
for n_model = 1:length(models)
    glm_specify_design_spmprep(sub_dir, models{n_model}, sprintf('testing_%02d', n_model), 1, 1)
end