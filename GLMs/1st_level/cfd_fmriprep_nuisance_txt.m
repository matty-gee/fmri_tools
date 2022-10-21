
function nuisance_txt = cfd_fmriprep_nuisance_txt(func_dir, nuisance)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creates a .txt from fmriprep confounds file with specified nuisance variables
    %
    % Arguments
    % ---------
    % func_dir : str
    %   path to find fmriprep confound file & save resulting text file
    % nuisance : str
    %   which nuisance variables to include separated by '+' signs
    %   OPTIONS: 'rp', 'rp24', 'fd', 'csf', 'csf4', 'acc_csf'
    %   example: 'rp+fd' --> realignment parameters + framewise displacement
    %   note: some confound variables don't have a value for 1st volume; use mean or 0
    %   
    % Returns
    % -------
    % str
    %   name of nuisance text file
    %
    % [By Matthew Schafer, github: @matty-gee; 2020ish] 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % minerva paths to add;  if not on minerva, no harm no foul
    addpath hpc/packages/minerva-centos7/spm/spm12
    addpath sc/arion/projects/k23/code/utilities
        
    f = filesep; % system specific
    split_ = split(func_dir, f);
    sub_id = split_{end-1};
    nuisance_txt = [func_dir f sub_id '_' nuisance '.txt'];
    
    if ~isfile(nuisance_txt)
        
        cfds_T     = struct2table(tdfread(spm_select('FPList', func_dir, '^s.*timeseries.tsv$')));
        cfd_struct = struct(); % add tables to a struct 
        
        %% motion parameters
        
        cfd_struct.rp   = table(cfds_T.trans_x, cfds_T.trans_y, cfds_T.trans_z, cfds_T.rot_x, cfds_T.rot_y, cfds_T.rot_z);
        
        cfd_struct.rp24 = table(cfds_T.trans_x, cfds_T.trans_y, cfds_T.trans_z,...% translations
                                 cfds_T.trans_x_power2, cfds_T.trans_y_power2, cfds_T.trans_z_power2,...% squared
                                 cfd_replace_firstvol(cfds_T.trans_x_derivative1, 'mean'),...% derivatives
                                 cfd_replace_firstvol(cfds_T.trans_y_derivative1, 'mean'),...
                                 cfd_replace_firstvol(cfds_T.trans_z_derivative1, 'mean'),...
                                 cfd_replace_firstvol(cfds_T.trans_x_derivative1_power2, 'mean'),...% derivatives squared
                                 cfd_replace_firstvol(cfds_T.trans_y_derivative1_power2, 'mean'),...
                                 cfd_replace_firstvol(cfds_T.trans_z_derivative1_power2, 'mean'),...
                                 cfds_T.rot_x, cfds_T.rot_y, cfds_T.rot_z,...% rotations
                                 cfds_T.rot_x_power2, cfds_T.rot_y_power2, cfds_T.rot_z_power2,...% squared 
                                 cfd_replace_firstvol(cfds_T.rot_x_derivative1, 'mean'),...% derivatives
                                 cfd_replace_firstvol(cfds_T.rot_y_derivative1, 'mean'),...
                                 cfd_replace_firstvol(cfds_T.rot_z_derivative1, 'mean'),...
                                 cfd_replace_firstvol(cfds_T.rot_x_derivative1_power2, 'mean'),...% derivatives squared
                                 cfd_replace_firstvol(cfds_T.rot_y_derivative1_power2, 'mean'),...
                                 cfd_replace_firstvol(cfds_T.rot_z_derivative1_power2, 'mean')); 
                         
        %% framewise displacement
        
        cfd_struct.fd = table(cfd_replace_firstvol(cfds_T.framewise_displacement, 'mean'));
        
        %% csf
        
        cfd_struct.csf  = table(cfds_T.csf);
                            
        cfd_struct.csf4 = table(cfds_T.csf,...
                                cfds_T.csf_power2,...
                                cfd_replace_firstvol(cfds_T.csf_derivative1, 'mean'),...
                                cfd_replace_firstvol(cfds_T.csf_derivative1_power2, 'mean'));
    
        %% aCompCor CSF PCs for up to 50% variance
        % wm, csf & wm+csf mask-based PCs that explain up to 50% variance are in confound file
        % -- fMRIPrep does high-pass filtering before running CompCor, so when using CompCor regressors include cosine_XX regressors
        % -- BUT if GLM estimation incls. separate hpf, do not include cosine_XX regressors in your design matrix...
        % -- Takeaway: if want to include acc in 1st level model, have to turn off hpf in SPM
        
        cfd_json  = cfd_json_read(spm_select('FPList', func_dir, '^s.*confounds_timeseries.json$'));
        cfd_names = fieldnames(cfds_T);
        acc_names = cfd_names(~cellfun(@isempty, regexp(cfd_names, '^a_comp_cor.*')));
        
        % only acompcor pcs estimated in CSF mask
        acc_csf_names = {}; 
        for a = 1 : length(acc_names)
            if strcmp(cfd_json.(acc_names{a}).Mask, 'CSF') 
                acc_csf_names{end + 1} = acc_names{a};
            end
        end
        acc_csf = cfds_T(:, acc_csf_names);

        % cosine regressors for hpf
        cfd_struct.acc_csf = [acc_csf, cfds_T(:, cfd_names(~cellfun(@isempty, regexp(cfd_names, '^cosine.*'))))];

        %% combine & output as txt file (maybe can use struct in spm?)
        
        nuisance_vars = strsplit(nuisance, '+');
        nuisance_T = [];
        for v = 1 : length(nuisance_vars)
            cfd_T = cfd_struct.(nuisance_vars{v});
            % have to change column names so dont overlap
            colnames = {};
            for c = 1 : width(cfd_T)
                colnames{c} = [nuisance_vars{v} num2str(c)]; 
            end
            cfd_T.Properties.VariableNames = colnames;
            nuisance_T = [nuisance_T, cfd_T];
        end
            
        writetable(nuisance_T, nuisance_txt, 'WriteVariableNames', 0);
 
    end 

end