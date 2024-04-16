%------------------------------------------------------------------------------------
% extract the parameters (betas, t-values, % signal change) from a masked functional image
% output in a .mat
% currently written for timeseries 
%
% NOTES: 
% - voxel dimensions should match - coregister(estimate & reslice) in spm first
% - use binarized masks
%------------------------------------------------------------------------------------

clear 

base_dir = '/Volumes/synapse/projects/SocialSpace/Projects';
out_dir  = [base_dir '/SNT-fmri_place/New_analyses/Timeseries'];

% subjects
preprc_dir = [base_dir '/SNT-fmri_CUD/Data/Scans/spm_preprocessing/subs'];
sub_dirs   = cellstr(spm_select('FPList', preprc_dir, 'dir', '^.*$'));
sub_dirs   = flip(sub_dirs, 1);

% masks
mask_dir = [base_dir '/SNT-fmri_place/Masks'];
masks = cellstr(spm_select('FPList', mask_dir, '^.*nii$'));

% filter mask list
include = {'.*HPC.*thr25.*', '.*M1.*thr25.*', '.*mPFC.*thr25.*'};
masks = masks(filter_cell_array_regex(masks, include));
exclude = {'.*bilateral.*', '.*-CA.*', '.*-DG.*', '.*-sub.*'};
masks = masks(~filter_cell_array_regex(masks, exclude));

%------------------------------------------------------------------------------------
% loop over subjects & load functional images
%------------------------------------------------------------------------------------

for s = 1 : length(sub_dirs)
    
    sub_dir = [sub_dirs{s} '/func'];
    [~, sub_id, ~] = fileparts(sub_dirs{s});

    out_fname = [out_dir '/' sub_id '_timeseries.mat'];
    if isfile(out_fname) 
        data = load(out_fname);  % load it if it exists
        data = data.data;
    else        
        data = struct(); 
    end

    % images could be '^wau.*nii$', '^spmT.*nii$', '^beta.*nii$'
    func_img_fnames = cellstr(spm_select('ExtFPList', sub_dir, '^wau.*nii$')); 

    %------------------------------------------------------------------------------------
    % loop over masks & extract timeseries

    % cross-platform compatibility test: 03.24.23 
    % - spm_mask_vol ensures nilearn created masks work same in matlab
        % matlab_mask  = [mask_dir '/Repl_L_mpfc_harvardoxford_maxprob-thr25-1mm_matlab.nii'];
        % nilearn_mask = [mask_dir '/Repl_L_mpfc_harvardoxford_maxprob-thr25-1mm.nii'];
        % isequal(load_mask(nilearn_mask), load_mask(matlab_mask))
    %------------------------------------------------------------------------------------

    for m = 1 : length(masks)

        % mask details
        mask_fname = masks{m};
        [~, mask_name, ~] = fileparts(mask_fname);
        split = strsplit(mask_name, '_');
        roi_name = [split{2} '_' split{3} '_' replace(replace(split{5}, 'maxprob-', ''), '-1mm', '')];
        
        % if roi_name exists in data, skip
        run = 1;
        for i = 1 : length(data)
           if ~isempty(fieldnames(data))
               if strcmp(data(i).roi_name, roi_name)
                   run = 0;
                   disp(['skipping ' sub_id ' ' roi_name])
               end
           end
        end
        
        if run
            
            disp(['running ' sub_id ' ' roi_name])
            
            % get mask indices
            mask_xyz = load_mask(mask_fname); % enforces compatibility between nilearn & spm

            % check that 3rd dimension isnt all same value
            if all(mask_xyz(3,:) ==  mask_xyz(3,1))
               disp('WARNING: 3rd dimension is all same value') 
            end

            % extract timeseries
            timeseries = get_mask_data(func_img_fnames, mask_xyz);

            % xyz coordinates
            xyz = array2table(mask_xyz');
            xyz.Properties.VariableNames = {'x','y','z'};

            % output in struct
            if isempty(fieldnames(data)), next = 1;
            else,                         next = length(data) + 1;
            end
            data(next).roi_name = roi_name;
            data(next).mask_fname = mask_fname;
            data(next).timeseries = timeseries;
            data(next).xyz = xyz;
            data(next).datetime = datestr(datetime);
            
        end 
        
    end

    save(out_fname, 'data')
end