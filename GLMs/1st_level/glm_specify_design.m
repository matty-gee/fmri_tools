function glm_specify_design(func_dir, model, glm_name, verbose)

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
% directories
%------------------------------------------------------------------------------


addpath /sc/arion/projects/k23/code/fmri_utilities % add circ_mean here
addpath /hpc/packages/minerva-centos7/spm/spm12

f = filesep; % system-specific 

% base_dir = [f 'sc' f 'arion' f 'projects' f 'k23'];
base_dir = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD';
split_   = strsplit(func_dir, '/');
sub_id   = erase(split_{end-1}, 'sub-P');

beh_dir  = [base_dir f 'Data' f 'SNT']; 
glm_dir  = [base_dir f 'GLMs' f glm_name f 'subs' f sub_id]; % for the output
if ~exist(glm_dir, 'dir'), mkdir(glm_dir), end
disp(['Preparing to run GLM for ' sub_id '...'])


%------------------------------------------------------------------------------
% images, masking & nuisance regrs
%------------------------------------------------------------------------------


% funcional images
func_name = 'wau_func'; 
% img_unzip(func_dir, func_fname)
func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_name '.nii$']));

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
% behavioral files to define onsets, durations, pmods
%------------------------------------------------------------------------------


behavior = readtable([beh_dir f 'Behavior' f 'snt_' sub_id '_behavior.xlsx']);
behavior = sortrows(behavior, 'decision_num', 'ascend');
timing   = readtable([beh_dir f 'Timing' f 'snt_' sub_id '_timing.xlsx']); % change to hardcode if single timing file
timing   = sortrows(timing, 'onset', 'ascend');


%------------------------------------------------------------------------------
% make & compute glm
%------------------------------------------------------------------------------   


% make desgin
glm     = glm_make_design(model, timing, behavior, glm_dir, verbose);
glm.tr  = 1; 
glm.mr  = 70; % microtime resolution: number of slices collected
glm.mo  = glm.mr/2; 
glm.hrf = [0 0]; 
glm.hpf = 128; % default = 128s (1/128 Hz) [or: max(diff(pmod_model.decision_onsets))]
glm.cvi = 'FAST'; % prewhitening to remove autocorrelated signal 

glm.model_post_decisions = 0;
glm.write_residuals      = 0; % write residuals (1) or not (0) [eg, 1 if want to do func conn after]

% clean up glm design object so can pass glm.cond directly into spm
glm.cond = rmfield(glm.cond, 'trials');
for n_cond = 1:length(glm.cond)
    glm.cond(n_cond).pmod = rmfield(glm.cond(n_cond).pmod, 'normalized');
end

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