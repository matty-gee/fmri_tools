%% Matt Schafer
% check motion, create a jpg w motion estimates

rp_dir = '/Users/matthew/Desktop/SocialSpace_Data/Scans/Raw/CUD/NEW';
cd(rp_dir);
files = dir('P*');
%files = dir('rp*.txt');
figure;

for f = 1:length(files)
    
    cd(rp_dir);
    cd([files(f).name '/Func']);
    subnum=files(f).name;
    rp=dir('rp*.txt');
    rp = spm_load(rp.name);
    motion_img = figure;
    subplot(2, 1, 1); plot(rp(:, 1:3));
    set(gca,'xlim',[0 size(rp, 1) + 1]);
    subplot(2, 1, 2); plot(rp(:, 4:6));
    set(gca,'xlim',[0 size(rp,1) + 1]);
    
    split = regexp(files(f).name, '[.]', 'split');
    filename = split{1};
    saveas(motion_img,filename,'jpg')
    
end


% spm_load(spm_select) % if want to select with gui