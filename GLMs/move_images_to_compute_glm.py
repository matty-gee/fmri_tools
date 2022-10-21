import glob, os, shutil

prprc_dir = '/sc/arion/projects/k23/derivatives_fieldmaps'
out_dir   = prprc_dir + '/for_glms'
if not os.path.exists(out_dir): os.makedirs(out_dir)

for sub_dir in glob.glob(prprc_dir + '/fmriprep/sub-P*'):
    if 'html' not in sub_dir:
        sub_id = sub_dir.split('/')[-1]
        print(sub_id)

        sub_out_dir = out_dir + '/' + sub_id
        if not os.path.exists(sub_out_dir): os.makedirs(sub_out_dir)

        func_img = glob.glob(sub_dir + '/func/*smoothed6.nii')[0]
        mask_img = glob.glob(sub_dir + '/func/*brain_mask.nii')[0]
        cfd_json = glob.glob(sub_dir + '/func/*confounds_timeseries.json')[0]
        cfd_tsv  = glob.glob(sub_dir + '/func/*confounds_timeseries.tsv')[0]

        for fname in [func_img, mask_img, cfd_json, cfd_tsv]:
            shutil.copy(fname, sub_out_dir + '/' + fname.split('/')[-1]) 

            # try:
            #     shutil.copy(fname, sub_out_dir + '/' + fname.split('/')[-1]) 
            # except:
            #     print('error w/ ' + fname)