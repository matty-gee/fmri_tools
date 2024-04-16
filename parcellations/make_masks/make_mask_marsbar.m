clear
out_dir = '/Users/matty_gee/Desktop/2D_place/Pmod_masks';

% define rois
radii = [5; 10; 15]; % in mm

roi_dir = '/Users/matty_gee/Dropbox/Projects/fmri_tools/parcellations';
roi_table = readtable([roi_dir '/social-space_rois_new.xlsx']);
% all_coords = [-45 41 -2; 51 50 1]; % coordinates are X,Y,Z
% region_names = {'L-IFG'; 'R-IFG'};

% make rois
fprintf('\n');
for r = 1 : length(radii)
    radius = radii(r); 
    for i = 1 : height(roi_table)

        x = roi_table.MNI_X(i);
        y = roi_table.MNI_Y(i);
        z = roi_table.MNI_Z(i);
        coords = [x,y,z];
        
        region_name = roi_table.ROI{i};
        study_name  = roi_table.Study{i};
        cond_name   = roi_table.Condition{i};

        fprintf('Working on ROI %d/%d...', i, height(roi_table));
        out_name = fullfile(out_dir, sprintf('%s_%s_%s_%i_%i_%i_%dmm', study_name, cond_name, region_name,...
                                                                       coords(1), coords(2), coords(3), radius));

        sphere_roi = maroi_sphere(struct('centre', coords, 'radius', radius)); % make sphere
    %     saveroi(sphere_roi, [out_name '.mat']); % save .mat
        save_as_image(sphere_roi, [out_name '.nii']); % save .nii

        fprintf('done.\n');
    end
end