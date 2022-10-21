clear 

%% data

home_dir     = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD';
glm_dir      = [home_dir '/Analyses/GLMs/Pmods/egocentric_angles/fmriprep-rp+fd+acc_csf/'];
subs_info    = readtable([home_dir '/participants_info.xlsx']);
summary_data = readtable([home_dir '/Data/Summary/All-data_summary_n78.xlsx']);
contrasts    = cellstr(spm_select('FPList', [glm_dir '/contrast_images'], '.*con_0003.nii'));
scans        = {};

%% covariates
% iCFI: interactions w/ experimental factor (creates an addtl regressor)
% - 1 = None
% iCC: centering (dont do for categorical
% - 1 = overall mean, 5 = no centering

covs = {'CocTox', 'CocCrave', 'CocPastMonth', 'CocFirstUse', 'CocYrsUse'};
cov = struct();
for c = 1 : length(covs)
    cov(c).cname = covs{c};
    cov(c).c     = [];
    cov(c).iCFI  = 1; % use flexible factorial if want group factor
    if strcmp(covs{c}, 'CocTox') || strcmp(covs{c}, 'Sex') || strcmp(covs{c}, 'White') % categoricals: dont mean center 
        cov(c).iCC = 0; 
    else
        cov(c).iCC = 1;
    end
end


%% contrast images etc
for s = 1 : length(contrasts)
    
    % get subj id from img file name
    split_ = strsplit(contrasts{s}, '/');
    split_ = strsplit(split_{end}, '_');
    sub_id = erase(split_{1}, 'sub-P');
    
    % check for inclusion
    sub_info = subs_info(any(ismember(subs_info.sub_id, str2double(sub_id)), 2), :);
    sub_data = summary_data(any(ismember(summary_data.sub_id, str2double(sub_id)), 2), :);    
    if (sub_info.incl == 1) && strcmp(sub_info.dx{1}, 'CD')
        
        scans{end+1, 1} = contrasts{s};
        
        % covariates
        for c = 1 : length(covs)
            
            switch covs{c}
               case 'Sex'
                    sex = sub_data.sex{1};
                    if strcmp(sex, 'female'), cov(c).c(end+1, 1) = 0;
                    else,                     cov(c).c(end+1, 1) = 1;
                    end
                case 'Age' 
                    cov(c).c(end+1, 1) = sub_data.age_years;
                case 'White' % binary
                    cov(c).c(end+1, 1) = sub_data.race___white; 
                case 'Edu' % ordinal
                    cov(c).c(end+1, 1) = str2num(sub_data.asi_education{1});  % coded weird in xlsx rn...
                case 'RT'
                    cov(c).c(end+1, 1) = sub_data.rt_mean;
                case 'Memory'
                    cov(c).c(end+1, 1) = sub_data.memory_mean;
                case 'MissedTrials'
                    cov(c).c(end+1, 1) = sub_data.missed_trials;
                case 'CocTox'  % toxicology on day of mri: objective state at scan
                    cov(c).c(end+1, 1) = sub_data.mri_utox___coc; 
                case 'CocCrave' % craving  on day of scan: subjective state at scan
                    cov(c).c(end+1, 1) = sub_data.cssa_4_coccrav_mri; 
                case 'CocPastMonth' % use recently (eg last month): longer time scale use
                    cov(c).c(end+1, 1) = sub_data.asi_coc_pastmonth; 
                case 'CocFirstUse' % age at first use
                    cov(c).c(end+1, 1) = sub_data.coc_age_1st_use;
                case 'CocYrsUse' % length of use
                    cov(c).c(end+1, 1) = sub_data.age_years - sub_data.coc_age_1st_use;
                case 'CTQ' % childhood trauma questionnaire
                    cov(c).c(end+1, 1) = sub_data.ctq_total_score_2;
                case 'RMSdiff' % root mean squared 
                    cov(c).c(end+1, 1) = sub_data.rms_diff_mean;
                case 'FD' % framewise displacement
                    cov(c).c(end+1, 1) = sub_data.fd;        
            end
        end
        
        % if cov is missing, remove the subj (prob cld b done at end, faster)
        for cc = 1 : length(covs)
            if isnan(cov(cc).c(end, 1))
                disp([sub_id ' has missing ' cov(cc).cname ' values - removing'])
                scans(end) = [];
                for ccc = 1 : length(covs), cov(ccc).c(end) = []; end
                continue 
            end
        end
    end
end

out_dir = [glm_dir '/angle_multiple-regression'];
if ~exist(out_dir, 'dir'), mkdir(out_dir), end
regression_compute(scans, cov, out_dir)