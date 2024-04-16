function save_V_as_nii(V, nii_name)
    % save V object from spm_vol as a nifti image, using marsbar
    save_as_image(maroi_image(V), nii_name)
