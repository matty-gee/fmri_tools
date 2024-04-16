% collect 3d images & concatenate into 4d
% for both beta & t-value imgs

clear

sub_list = cellstr(spm_select('List', cd, 'dir'));

for s = 1 : length(sub_list)

    cd(sub_list{s})
    beta_imgs_3d = spm_select('FPList', cd, '^beta.*nii$');
    spm_file_merge(beta_imgs_3d(1:63,:), 'merged_beta_4d.nii');

    tval_imgs_3d = spm_select('FPList', cd, '^spmT.*nii$');
    spm_file_merge(tval_imgs_3d(1:63,:), 'merged_tval_4d.nii');

    cd .. 

end


