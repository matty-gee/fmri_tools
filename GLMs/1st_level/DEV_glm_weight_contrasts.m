clear 

addpath /sc/arion/projects/k23/code/matlab_utilities
addpath /hpc/packages/minerva-centos7/spm/spm12

% TODO: maybe load in the SPM.mat or GLM.mat to read in the names?

f = filesep;

glm_dir = '/Users/matty_gee/Desktop/2D_place/Samples/Validation/GLMs/angle/subs';
d = dir(glm_dir); % get the folder contents
sub_dirs = d([d(:).isdir]); % remove all files (isdir property is 0)
sub_dirs = sub_dirs(~ismember({sub_dirs(:).name},{'.','..'})); % remove '.' and '..' 

for s = 1 : length(sub_dirs)
    
    sub_dir = [sub_dirs(s).folder '/' sub_dirs(s).name];

    % check if the contrasts are weighted 
    if ~isfile([sub_dir f 'con_0001.nii'])

        %------------------------------------------------------------------------------------
        % specify weights
        %------------------------------------------------------------------------------------

        glm_design = struct();

        % weight individual regrs (1) against baseline
        regrs = {'narrative'; 'decisions'; 'decisions*pov3d_angle'};
        n_regrs = length(regrs);
        tcons = {};
        for n_regr = 1:n_regrs

            tcon_weights = zeros(1, length(regrs));
            tcon_weights(n_regr) = 1;

            % positive 
            tcon_name = [regrs{n_regr} '+'];
            tcons = [tcons, tcon_name];
            glm_design.consess{length(tcons)}.tcon = struct('name', tcon_name, 'weights',  tcon_weights);

            % negative
            tcon_name = [regrs{n_regr} '-'];
            tcons = [tcons, tcon_name];
            glm_design.consess{length(tcons)}.tcon = struct('name', tcon_name, 'weights', -tcon_weights);

        end
        n_tcons = length(tcons);

        %------------------------------------------------------------------------------------
        % weight contrasts
        %------------------------------------------------------------------------------------


        batch{1}.spm.stats.con.spmmat                 = {[sub_dir f 'SPM.mat']};
        batch{1}.spm.stats.con.consess                = glm_design.consess;


        %------------------------------------------------------------------------------------
        % compute
        %------------------------------------------------------------------------------------


        try
            spm_jobman('run', batch);
        catch
            disp('Problem running spm_jobman. Running in interactive mode.')
            spm_jobman('interactive', batch);
        end
    end
end