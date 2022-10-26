function fmriprep_cfd_wmcsf_regressors(sub_dir, preprocessing)
    % OLD: dont use....
    % change this to accept arg for which func image to get wm & csf from
        
    addpath /hpc/packages/minerva-centos7/spm/spm12
    
    anat_dir = [sub_dir '/anat/'];
    func_dir = [sub_dir '/func/'];

    % load func images
    if strcmp(preprocessing, 'ica_aroma')
        func_fpattern = 'AROMAnonaggr_bold.nii';
    elseif strcmp(preprocessing, 'fmriprep')
        func_fpattern = 'preproc_bold.nii';
    end
    if ~isfile(spm_select('FPList', func_dir, ['^.*' func_fpattern '$']))
        gunzip(spm_select('FPList', func_dir, ['^.*' func_fpattern '.gz$']))
    end
    func_nii = cellstr(spm_select('ExtFPList', func_dir, ['^.*' func_fpattern '$']));
    
    %% make masks: unzip & binarize probability imgs at 95% & resample to functional images
    
    tpm_fnames = {'^sub.*MNI.*WM_probseg'; '^sub.*MNI.*CSF_probseg'};
    for n = 1 : length(tpm_fnames)

        % 1. unzip
        unzip_nifti(anat_dir, [tpm_fnames{n} '.nii$'])
%         if ~isfile(spm_select('FPList', anat_dir, [tpm_fnames{n} '.nii$']))           
%             gunzip(spm_select('FPList', anat_dir, [tpm_fnames{n} '.nii.gz$']))
%         end

        % 2. binarize anatomical
        bin_nii = [tpm_fnames{n} '_binarized-thresh95.nii$'];
        if ~isfile(spm_select('FPList', anat_dir, bin_nii))
            disp('Binarizing images')
            binarize_img(spm_select('FPList', anat_dir, [tpm_fnames{n} '.nii$']), 0.95)
        else
            disp('Images already binarized')
        end

        % 3. resample to functional image
        resamp_nii = [tpm_fnames{n} '_binarized-thresh95_resampled-func.nii$'];
        if ~isfile(spm_select('FPList', anat_dir, resamp_nii))
            disp('Resampling images')
            resample_img(spm_select('FPList', anat_dir, bin_nii), func_nii{1}, 'resampled-func');  
         else
            disp('Images already resampled to functional')
        end

    end

    %% apply resulting masks to the nonaggr bold nifti & get mean value per volume
    % -- dont have to rebinarize in spm

    [~, func_fname, ~] = fileparts(func_nii{1});
    if ~isfile(spm_select('FPList', func_dir, [func_fname '_wmcsf-thr95.csv']))

        disp('Computing mean wm & csf signals for functional')
        [csf_V, ~] = spm_read_vols(spm_vol(spm_select('FPList', anat_dir, '^sub.*MNI.*CSF.*resampled-func.nii$')), 1);
        [wm_V, ~]  = spm_read_vols(spm_vol(spm_select('FPList', anat_dir, '^sub.*MNI.*WM.*resampled-func.nii$')), 1);    

        % get mean per volume
        wm_95 = []; csf_95 = []; 
        for i = 1 : length(func_nii)
            v = spm_read_vols(spm_vol(func_nii{i})); % load in volume
            wm_95  = [wm_95; spm_global(wm_V .* v)];
            csf_95 = [csf_95; spm_global(csf_V .* v)];
        end

        % output csv & txt
        nuisance_table = table(wm_95, csf_95);
        writetable(nuisance_table, [func_dir '/' func_fname '_wmcsf-thr95.csv']);
        writetable(nuisance_table, [func_dir '/' func_fname '_wmcsf-thr95.txt'], 'WriteVariableNames', 0)
    end
 
end
