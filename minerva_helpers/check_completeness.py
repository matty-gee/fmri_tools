# check which of preprocessed subs dont have some folder in some analysis...
import glob, os

base_dir = '/sc/arion/projects/OlfMem/mgs/2D_place'
for sample in ['Validation', 'Initial']:
    preprc_dir = f'{base_dir}/Samples/{sample}/spm_preprocessing/subs'
    analysis_dir = f'{base_dir}/Samples/{sample}/GLMs/lsa_decision/subs'
    sub_dirs = [sub for sub in glob.glob(f'{preprc_dir}/*') if 'html' not in sub] # exclude html files
    for sub_dir in sub_dirs: 
        sub_id = sub_dir.split('/')[-1]
        if not os.path.exists(f'{analysis_dir}/{sub_id}/beta_0001.nii'):
            print(f'{sub_id} is missing from {analysis_dir}')