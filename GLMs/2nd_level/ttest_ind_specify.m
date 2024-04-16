clear 

% home_dir = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD';
home_dir = '/Users/matty_gee/Desktop/CUD';
glm_dir = [home_dir '/GLMs/Pmods/egocentric_angles/fmriprep-rp+fd/'];

% organize the contrast files...
subs_info = readtable([home_dir '/participants_info.xlsx']);
for cn = 1 : 2

    contrasts = cellstr(spm_select('FPList', [glm_dir '/contrast_images'], ['.*con_000' num2str(cn) '.nii']));
    cd_scans = {}; hc_scans = {};
    for ci = 1 : length(contrasts)

        % check for inclusion
        split_ = strsplit(contrasts{ci}, '/');
        split_ = strsplit(split_{end}, '_');
        sub_id = erase(split_{1}, 'sub-P');
        sub_info = subs_info(any(ismember(subs_info.sub_id, str2double(sub_id)), 2), :);
        if sub_info.incl == 1 
            if strcmp(sub_info.dx{1}, 'CD')
                cd_scans{end+1, 1} = contrasts{ci};
            else
                hc_scans{end+1, 1} = contrasts{ci};
            end
        end
    end

    if cn == 1
        out_dir = [glm_dir '/decision_ttest'];
    else
        out_dir = [glm_dir '/angle_ttest'];
    end

    if ~exist(out_dir, 'dir'), mkdir(out_dir), end
    ttest_ind(hc_scans, cd_scans, out_dir)
        
end