function xyz = load_mask(mask_fname)
    % mask_fname should be a string
    V         = spm_mask_vol(mask_fname); % returns data structure
    vols      = spm_read_vols(V, 1); % read in binary mask
    [x, y, z] = ind2sub(size(vols), find(vols>0)); % turn 1s into xyz coordinates
    xyz       = [x y z]'; % transpose coords: 3 dimensions by n voxels index matrix
    % warn if xyz is empty...?
end
