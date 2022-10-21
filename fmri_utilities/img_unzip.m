function img_unzip(img_dir, img_stem)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unzip an image, if needed
    %
    % Arguments
    % ---------
    % img_dir : str
    % img_stem : str
    %   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    addpath /hpc/packages/minerva-centos7/spm/spm12
    
    if isempty(spm_select('FPList', img_dir, [img_stem '.nii$']))
        disp('Unzipping image')
        gunzip(spm_select('FPList', img_dir, [img_stem '.nii.gz$'])); 
    else
        disp('Already unzipped')
    end    
end