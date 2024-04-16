function V = spm_mask_vol(mask_fname)

    % read in mask volume with spm_vol & ensure cross-platform (nilearn & spm) compatibility
    % - e.g., wrt datatype & nans
    
    V = spm_vol(mask_fname); % info in a struct
    
    % plane info
    V.pinfo = [1; 0; 352]; 
    
    % datatype
    info = niftiinfo(mask_fname);
    switch info.Datatype
        case 'uint8'
            d = 2;
        case 'int16'
            d = 4;
        case 'int32'
            d = 8;
        case {'float32', 'single'}
            d = 16;
        case {'int64', 'float64', 'double'}
            d = 64;
        case 'int8'
            d = 256;
        case 'uint16'
            d = 512;
        case 'uint32'
            d = 768;
    end
    V.dt = [d 0];
    

    % if need to save this as an image
    % save_as_image(maroi_image(V_out), 'mask_test.nii')
    
end