function data = get_mask_data(imgs, mask_xyz)
    
    %----------------------------------------------------------------
    % return the data from a mask
    % mask_xyz should be indices in 3 dimensions by n voxels
    %----------------------------------------------------------------
    
    data = zeros(length(imgs), length(mask_xyz)); % pre-allocate
    for i = 1 : length(imgs)  
        data(i, :) = spm_get_data(imgs{i}, mask_xyz)';      
    end
    
    % checks
    if length(mask_xyz) ~= size(data, 2)
        error('The voxel counts do not match')
    end     