% art repair requires 3d images
% normally would do this with smoothed images right before GLM estimation
% but this is for mvpa with unsmoothed images

clear 
if ismac
   base_dir = '/volumes/synapse/projects/SocialSpace';
elseif isunix 
   base_dir = '/mnt/synapse/projects/SocialSpace';
end
scans_dir = [base_dir '/Data/Scans/RT/fmriprep/derivatives/fmriprep']; 
scan_dirs = cellstr(spm_select('FPList', scans_dir, 'dir')); 

for s = 1 : length(scan_dirs)
    
    %% directories
    
    sub_dir = scan_dirs{s};
    pathparts = strsplit(sub_dir, '/');
    sub_id = pathparts{end};
    preprc_dir = [sub_dir '/func'];
    artrepr_dir = [preprc_dir '/art_repair_unsmoothed'];    
    subparts = strsplit(sub_id, '_');
        
    %% check if the repaired file exists

    if isfile([artrepr_dir '/art_repaired.txt'])
        disp([sub_id ' is repaired. Going to next subject.']) 

    else

        cd(sub_dir)
        if ~exist(artrepr_dir, 'dir'), mkdir(artrepr_dir); end

        %% 4D -> 3D

        % (normally would be with smoothed files but not using smoothed for rsa)
        spm_file_split(spm_select('FPList', preprc_dir, '^w.*nii$'), artrepr_dir);
        unrepaired_3d = spm_select('FPList', artrepr_dir, '^w.*nii$');

        %% root mean square realignment estimates
        if isempty(spm_select('FPList', artrepr_dir, 'art_rms.txt'))
            % RMS percentage, actual RMS, image mean over mask
            var_stats = art_rms(unrepaired_3d);
            rms_T = array2table(var_stats);
            rms_T.Properties.VariableNames = {'rms_percentage', 'rms', 'image_mean'};
            writetable(rms_T, [artrepr_dir '/art_rms.txt']);
        end

        %% find & deweight bad volumes

        % FORMAT art_global(Images, RealignmentFile, HeadMaskType, RepairType)
        %    Images  = Full path name of images to be repaired.
        %       For multiple sessions, all Images are in one array, e.g.SPM.xY.P
        %    RealignmentFile = Full path name of realignment file
        %       For multiple sessions, cell array with one name per session.
        %    HeadMaskType  = 1 for SPM mask, = 4 for Automask
        %    RepairType = 1 for ArtifactRepair alone (0.5 movement and add margin).
        %               = 2 for Movement Adjusted images  (0.5 movement, no margin)
        %               = 0 No repairs are done, bad scans are found.
        %                   Listed in art_suspects.txt for motion adjustment.
        
        % create rp file from fmriprep confound file
        confounds = tdfread(spm_select('FPList', preprc_dir, '^s.*confounds.*.tsv$'));
        writetable(table(confounds.trans_x, confounds.trans_y, confounds.trans_z, confounds.rot_x, confounds.rot_y, confounds.rot_z), [artrepr_dir '/rp.txt'], 'WriteVariableNames', 0);

        % art repair stuff 
        if isempty(spm_select('FPList', artrepr_dir, 'art_repaired.txt'))
            
            rp_file = spm_select('FPList', artrepr_dir, 'rp.txt');
            head_mask_type = 4;
            repair_type = 2;
            art_global(unrepaired_3d, rp_file, head_mask_type, repair_type)

            %% clean up

            repaired_3d = spm_select('FPList', artrepr_dir, '^v.*nii$');
            [~, fname, ~] = fileparts(repaired_3d(1,1:end));
            fname_split = strsplit(fname, '_');
            spm_file_merge(repaired_3d, [artrepr_dir '/' fname_split{1} '_' fname_split{2} '.nii']);
            for i=1:length(repaired_3d(:,1)), delete(repaired_3d(i,1:end)), end
            for i=1:length(unrepaired_3d(:,1)), delete(unrepaired_3d(i,1:end)), end
        end
    end
end
