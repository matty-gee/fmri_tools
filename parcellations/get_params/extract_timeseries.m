function extract_timeseries(sub_dir)

    %------------------------------------------------------------------------------------
    % loop over masks & extract timeseries
    % doesnt make any assumptions about data (e.g. doesnt need to know tr etc)

    % cross-platform compatibility test: 03.24.23 
    % - spm_mask_vol ensures nilearn created masks work same in matlab
        % matlab_mask  = [mask_dir '/Repl_L_mpfc_harvardoxford_maxprob-thr25-1mm_matlab.nii'];
        % nilearn_mask = [mask_dir '/Repl_L_mpfc_harvardoxford_maxprob-thr25-1mm.nii'];
        % isequal(load_mask(nilearn_mask), load_mask(matlab_mask))
    %------------------------------------------------------------------------------------
    
    addpath /sc/arion/projects/OlfMem/mgs/2D_place/Code/matlab_utilities
    addpath /sc/arion/projects/OlfMem/mgs/2D_place/Code/extract_timeseries
    addpath /hpc/packages/minerva-centos7/spm/spm12

    [~, sub_id, ~] = fileparts(sub_dir);
    disp(['Extracting timeseries for ' sub_id])
    
    base_dir = '/sc/arion/projects/OlfMem/mgs/2D_place';
    out_dir  = [base_dir '/Trajectory_analyses/timeseries/subs'];

    % load the timeseries file, if it exists
    out_fname = [out_dir '/' sub_id '_timeseries.mat'];
    if isfile(out_fname) 
        data = load(out_fname); % load it if it exists
        data = data.data;
    else        
        data = struct(); 
    end

    % images could be '^wau.*nii$', '^spmT.*nii$', '^beta.*nii$'
    func_img_fnames = cellstr(spm_select('ExtFPList', [sub_dir '/func'], 'wau_func.nii'));
    fprintf('Found %d images\n', length(func_img_fnames))

    % masks
    mask_dir = [base_dir '/Masks'];
    masks = cellstr(spm_select('FPList', mask_dir, '^.*nii$'));
    % exclude = {'.*bilateral.*'};
    % masks = masks(~filter_cell_array_regex(masks, exclude));
    fprintf('Found %d masks\n',length(masks))

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
                   disp(['already extracted ' roi_name])
               end
           end
        end

        if run

            disp(['extracting ' roi_name])

            try 
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
                
            catch
                
                disp(['error: skipping ' roi_name])
                
            end
        end 

    end

    save(out_fname, 'data')