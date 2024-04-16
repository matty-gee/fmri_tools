function vbm_run(sub_dir)

    addpath /hpc/packages/minerva-centos7/spm/spm12
    addpath /sc/arion/projects/k23/code/utilities

    gm_image = [sub_dir '/mwp1' sub_dirs(s).name '_T1w.nii'];
    if ~isfile(gm_image)
        
        % segment
        t1_image = spm_select('FPList', sub_dir, '.*T1w.nii.gz');
        cat_segmentation(t1_image) % will also decompress images
        
        % smooth
        img_smooth(gm_image, 8) % smoothing at 8mm fwhm
    end  