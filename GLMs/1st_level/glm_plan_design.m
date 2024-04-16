function glm_plan_design(sub_dir, preprc, model, hpf, write_residuals, verbose, debug)

%------------------------------------------------------------------------------
% Specify the design for a first-level GLM
% 
%
% Arguments
% ---------
% sub_dir: overall directory for subjects's functional images
%
% [By Matthew Schafer, github: @matty-gee; 2020ish] 
%------------------------------------------------------------------------------

% TODO: think about the architecture here
% add an LSA all option that will estimate every single trial

%------------------------------------------------------------------------------
% directories
%------------------------------------------------------------------------------


addpath /sc/arion/projects/OlfMem/mgs/2D_place/Code/matlab_utilities
addpath /sc/arion/projects/OlfMem/mgs/2D_place/Code
addpath /hpc/packages/minerva-centos7/spm/spm12

f = filesep; % system-specific 

func_dir = [sub_dir f 'func'];
split    = strsplit(sub_dir, '/');
sub_id   = erase(split{end}, 'sub-P');

if debug % debug locally
    
    beh_dir  = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/Data/SNT'; 
    glm_dir  = ['/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD' f 'TEST_GLM'];
    
elseif contains(pwd, '/Users/matthew') % local desktop
    
    base_dir = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD';
    
elseif contains(pwd, '/Users/matty_gee') % local laptop
    
    base_dir = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD';
    beh_dir  = [base_dir f 'Data' f 'SNT']; 
    glm_dir  = [base_dir f 'Analyses' f 'GLMs' f preprc '_' model f 'subs' f sub_id];
    
elseif contains(pwd, '/sc/arion/projects') % minerva

    base_dir = [f 'sc' f 'arion' f 'projects' f 'OlfMem' f 'mgs' f '2D_place'];
    sample   = split{end-3};
    beh_dir  = [base_dir f 'Samples' f sample f 'SNT']; 
    glm_dir  = [base_dir f 'Samples' f sample f 'GLMs' f model '_' num2str(hpf) 'hpf' f 'subs' f sub_id];
    
end

if ~exist(glm_dir, 'dir'), mkdir(glm_dir), end
disp(['Preparing to run GLM for ' sub_id '...'])


%------------------------------------------------------------------------------
% behavioral files to define onsets, durations, pmods
%------------------------------------------------------------------------------


behavior = readtable([beh_dir f 'Behavior' f 'SNT_' sub_id '_behavior.xlsx']);
behavior = sortrows(behavior, 'decision_num', 'ascend');
timing   = readtable([beh_dir f 'Timing' f 'SNT_' sub_id '_timing.xlsx']); % change to hardcode if single timing file
timing   = sortrows(timing, 'onset', 'ascend');

if verbose
    fprintf('%d timing file size\n', size(timing))
    fprintf('%d behavioral file size\n', size(behavior))
end


%------------------------------------------------------------------------------
% images, masking & nuisance regrs
%------------------------------------------------------------------------------


% funcional images
if strcmp(preprc, 'spm')
    func_fname = 'wau_func';  % func_name = '^wau.*';
    func_imgs  = cellstr(spm_select('ExtFPList', func_dir, [func_fname '.nii']));
elseif strcmp(preprc, 'fmriprep')
    func_fname = '^s.*preproc_bold';
    img_unzip(spm_select('FPList', func_dir, [func_fname '.nii$']))
    func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_fname '.nii$']));
end    

% if not lsa, smooth 
if (~strcmp(model, 'lsa_onset')) && (~strcmp(model, 'lsa_decision')) && (~debug) 
    img_smooth(spm_select('ExtFPList', func_dir, [func_fname '.nii']), 6)
    func_imgs = cellstr(spm_select('ExtFPList', func_dir, [func_fname '_smoothed6.nii']));  
end    

% nuisance regrs
if strcmp(preprc, 'spm')
    nuisance_txt = spm_select('FPList', func_dir, 'rp.txt');
elseif strcmp(preprc, 'fmriprep')
    % TODO: need to define cfds
    nuisance_txt = fmriprep_cfd_nuisance_txt(func_dir, cfds);
end 


%------------------------------------------------------------------------------
% make design & compute glm
%------------------------------------------------------------------------------   


glm = glm_make_design(model, timing, behavior, glm_dir, verbose);

% estimating the events [detect automatically b/c initial is diff. from validation]
if length(func_imgs) == 1570
    glm.tr  = 1; 
    glm.mr  = 70; % microtime resolution: number of slices collected
else
    glm.tr  = 2; 
    glm.mr  = 36;
end
glm.mo  = glm.mr/2; 
glm.hrf = [0 0]; 

% filtering etc
glm.hpf = hpf; % high-pass filter in seconds; default = 128s (1/128 Hz)
glm.cvi = 'FAST'; % prewhitening to remove autocorrelated signal 

% masking
if strcmp(preprc, 'spm')
    glm.mthresh  = 0.50; % default=0.8
    glm.mask_img = '';
elseif strcmp(preprc, 'fmriprep')
    img_unzip(func_dir, '^s.*brain_mask')  
    glm.mthresh  = -Inf; % -Inf allows for explicit masks
    glm.mask_img = spm_select('FPList', func_dir, '^s.*brain_mask.nii$');
end

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