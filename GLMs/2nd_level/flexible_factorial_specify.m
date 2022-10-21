clear 

% ADD DOCSTRINGS...

% specify design matrix flexibly
% -- 3 stages: specify factors, add scans, build design matri
% -- best for anovas where dont want to test all possible main & interaction fx

%% data

home_dir     = '/Volumes/synapse/projects/SocialSpace/Projects/SNT-fmri_CUD';
glm_dir      = [home_dir '/Analyses/GLMs_fieldmaps_rp/angle_decision/'];

subs_info    = dir([home_dir '/participants_qc_n*']);
subs_info    = readtable([home_dir '/' subs_info.name]);
summary_data = dir([home_dir '/Data/Summary/All-data_summary_n*']);
summary_data = readtable([home_dir '/Data/Summary/' summary_data.name]);

contrasts    = cellstr(spm_select('FPList', [glm_dir '/images'], '.*con_0003.nii')); 

subjects     = struct();
out_dir      = [glm_dir '/flexible-factorial'];
if ~exist(out_dir, 'dir'), mkdir(out_dir), end

%% covariates
% iCFI: interactions w/ experimental factor (creates an addtl regressor)
% - 1 = None, 2-4 = w/ factor 1-3
% iCC: centering
% - 1 = overall mean, 2-4 = factor 1-3 mean, 5 = no centering

covs = {'CTQ', 'Age', 'Sex'}; 
cov  = struct();
for c = 1 : length(covs)
    
    cov(c).cname = covs{c};
    cov(c).c = [];
    
    % interaction w/ group factor: creates a predictor w/ each level of factor
    if strcmp(covs{c}, 'CTQ')
        cov(c).iCFI = 2; % interaction
    else
        cov(c).iCFI = 1;
    end
    
    % mean center 
    if any(strcmp(covs{c}, {'CocTox','Sex','White'})) % categoricals: dont 
        cov(c).iCC = 0; 
    else
        cov(c).iCC = 1; % should i mean center w/in group for the others?
    end
    
end

%% organize 

n = 0; n_cd = 0; n_hc = 0;
for s = 1 : length(contrasts)
    
    % check for inclusion
    split_   = strsplit(contrasts{s}, '/');
    split_   = strsplit(split_{end}, '_');
    sub_id   = erase(split_{1}, 'sub-P');
    sub_info = subs_info(any(ismember(subs_info.sub_id, str2double(sub_id)), 2), :);
    sub_data = summary_data(any(ismember(summary_data.sub_id, str2double(sub_id)), 2), :);
    
    % if (sub_info.memory_incl == 1) && (sub_info.other_incl == 1) 
    if sub_info.memory_incl == 1
        
        % get subject & define their group 
        subjects(n+1, 1).scans = contrasts{s};
        if strcmp(sub_info.dx{1}, 'HC')
            subjects(n+1, 1).conds = 1; % group factor
            n_hc = n_hc + 1;
        else
            subjects(n+1, 1).conds = 2;
            n_cd = n_cd + 1;
        end
        
        % get covariates
        for c = 1 : length(covs)

            switch covs{c}
               case 'Sex'
                    sex = sub_data.sex{1};
                    if strcmp(sex, 'female'), cov(c).c(n+1, 1) = 0;
                    else,                     cov(c).c(n+1, 1) = 1;
                    end
                case 'Age' 
                    cov(c).c(n+1, 1) = sub_data.age_years;
                case 'White' % binary
                    cov(c).c(n+1, 1) = sub_data.race___white; 
                case 'Edu' % ordinal
                    cov(c).c(n+1, 1) = str2num(sub_data.asi_education{1});  % coded weird in xlsx rn...
                case 'RT'
                    cov(c).c(n+1, 1) = sub_data.rt_mean;
                case 'Memory'
                    cov(c).c(n+1, 1) = sub_data.memory_mean;
                case 'MissedTrials'
                    cov(c).c(n+1, 1) = sub_data.missed_trials;
                case 'CocTox'  % toxicology on day of mri: objective state at scan
                    cov(c).c(n+1, 1) = sub_data.mri_utox___coc; 
                case 'CocCrave' % craving  on day of scan: subjective state at scan
                    cov(c).c(n+1, 1) = sub_data.cssa_4_coccrav_mri; 
                case 'CocPastMonth' % use recently (eg last month): longer time scale use
                    cov(c).c(n+1, 1) = sub_data.asi_coc_pastmonth; 
                case 'CocFirstUse' % age at first use
                    cov(c).c(n+1, 1) = sub_data.coc_age_1st_use; 
                case 'CTQ' % childhood trauma questionnaire --> might want to turn into a factor?
                    cov(c).c(n+1, 1) = sub_data.ctq_total_score_2;
                case 'FD' % framewise displacement
                    cov(c).c(n+1, 1) = sub_data.fd_mean;        
            end
        end
        
        % if cov is missing, remove the subj (prob cld b done at end, faster)
        for cc = 1 : length(covs)
            if isnan(cov(cc).c(n+1, 1))
                disp([sub_id ' has missing ' cov(cc).cname ' values - removing'])
                subjects(n+1) = [];
                for ccc = 1 : length(covs), cov(ccc).c(n+1) = []; end
                continue 
            end 
        end 
        n = n + 1; % if sub is included
    end
end

% turn CTQ into categorical
if any(strcmp(covs, 'CTQ'))
    disp('Converting CTQ into categorical with median split')
    i = find(strcmp(covs, 'CTQ')==1);
    cov(i).c = (cov(i).c > median(cov(i).c)) * 1; % code > median as 1
end

%% run

disp(['HC n=' num2str(n_hc) ', CD n=' num2str(n_cd)])
flexible_factorial_compute(subjects, cov, out_dir)
cd(out_dir)