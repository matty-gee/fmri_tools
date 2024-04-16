function run_glm_jobs()

    addpath /sc/arion/projects/OlfMem/mgs/2D_place/Code/matlab_utilities
    addpath /sc/arion/projects/OlfMem/mgs/2D_place/Code
    addpath /hpc/packages/minerva-centos7/spm/spm12

    base_dir =  '/sc/arion/projects/OlfMem/mgs/2D_place/Samples';
    samples  = {'Initial', 'Validation'};
    models   = {'lsa_onset', 'lsa_decision', 'dimensions', 'cos_angle', 'distance', 'sin_angle', 'polar', 'character'};
    hpfs     = {128, 250};
    
    for mod = 1 : length(models)
        model = models{mod};
        for hp = 1 : length(hpfs)
            hpf = hpfs{hp};
            for samp = 1 : length(samples)
                sample = samples{samp};

                % find all sub_dirs in preprc_dir
                sub_dirs = cellstr(spm_select('FPList', [base_dir '/' sample '/spm_preprocessing/subs'], 'dir'));
                disp(['Running ' num2str(length(sub_dirs)) ' for ' model])

                % loop over
                for s = 1 : length(sub_dirs)
                    disp(['Running ' sub_dirs{s}])
                    glm_plan_design(sub_dirs{s}, 'spm', model, hpf, 0, 0, 0)
                end
                
            end
        end
    end
    