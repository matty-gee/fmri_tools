%% Matt Schafer
% create a jpg w motion parameter estimates

clear

rp_dir  = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/QC/motion_parameters/txts';
img_dir = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD/QC/motion_parameters/imgs';
files = dir(fullfile(cd, '**', '*rp.txt'));

for f = 1:length(files)
    
    rp = spm_load([files(f).folder '/' files(f).name]);
    
    motion_img = figure('Visible', 'off');
    
    % xyz translations
    subplot(2, 2, 1);
    plot(rp(:, 1:3));
    set(gca,'xlim', [0 size(rp, 1) + 1]);
    set(gca,'ylim', [-4 4]); % constant range to compare subjects
    title('XYZ translations [-4, 4]')

    subplot(2, 2, 2);
    plot(rp(:, 1:3));
    set(gca,'xlim', [0 size(rp, 1) + 1]);
    title('XYZ translations')
    
    
    % xyz rotations
    subplot(2, 2, 3); 
    plot(rp(:, 4:6));
    set(gca,'xlim', [0 size(rp, 1) + 1]);
    set(gca,'ylim', [-0.1, 0.1]);
    title('XYZ rotations [-0.1, 0.1]')
  
    subplot(2, 2, 4); 
    plot(rp(:, 4:6));
    set(gca,'xlim', [0 size(rp, 1) + 1]);
    title('XYZ rotations')
    
    split = regexp(files(f).name, '[.]', 'split');
    saveas(motion_img, [img_dir '\' split{1}], 'jpg')
    
end