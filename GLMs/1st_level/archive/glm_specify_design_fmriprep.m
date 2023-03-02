function glm_specify_design_fmriprep(func_dir, model, glm_name, write_residuals, verbose, debug)

%------------------------------------------------------------------------------
% Specify the design for a first-level GLM
% 
%
% Arguments
% ---------
% func_dir : str
% model : str or cell array (see below)
% glm_name : str
% write_residuals : bool
% verbose : bool
% 
%
% [By Matthew Schafer, github: @matty-gee; 2020ish] 
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
% directories
%------------------------------------------------------------------------------


addpath /sc/arion/projects/k23/code/matlab_utilities 
addpath /sc/arion/projects/k23/code/GLMs
addpath /hpc/packages/minerva-centos7/spm/spm12

f = filesep; % system-specific 

base_dir = [f 'sc' f 'arion' f 'projects' f 'k23'];
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
func_fname = '^s.*preproc_bold';
img_unzip(func_dir, func_fname)
func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_name '.nii$']));

if strcmp(glm_name, 'lsa') 
    func_imgs = cellstr(spm_select('ExtFPList', func_dir, '^s.*preproc_bold.nii$')); % use unsmoothed images
else
    img_smooth(spm_select('ExtFPList', func_dir, [func_fname '.nii$']), 6)
    func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_fname '_smoothed6.nii$']));  
end

% nuisance regrs
nuisance_txt = fmriprep_cfd_nuisance_txt(func_dir, cfds);


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


glm = glm_make_design(model, timing, behavior, glm_dir, verbose);

% estimating the events 
glm.tr  = 1; 
glm.mr  = 70; % microtime resolution: number of slices collected
glm.mo  = glm.mr/2; 
glm.hrf = [0 0]; 

% filtering etc
glm.hpf = 128; % default = 128s (1/128 Hz) [or: max(diff(pmod_model.decision_onsets))]
glm.cvi = 'FAST'; % prewhitening to remove autocorrelated signal 

% masking
img_unzip(func_dir, '^s.*brain_mask')  
glm.mthresh  = -Inf; % -Inf allows for explicit masks
glm.mask_img = spm_select('FPList', func_dir, '^s.*brain_mask.nii$');

% write out residuals (e.g. for fc)
glm.write_residuals = write_residuals; 

% clean up glm design object so can pass glm.cond directly into spm
glm.cond = rmfield(glm.cond, 'trials');
for n_cond = 1:length(glm.cond)
    glm.cond(n_cond).pmod = rmfield(glm.cond(n_cond).pmod, 'normalized');
end

if ~debug, glm_compute(func_imgs, glm, nuisance_txt, glm_dir), end


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