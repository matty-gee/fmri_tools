clear 

base_dir = '/Users/matty_gee/Desktop/2D_place/Samples';
img_dirs = {[base_dir '/Initial/GLMs/polar/contrast_images'], 
            [base_dir '/Validation/GLMs/polar/contrast_images']};

for d = 1 : length(img_dirs)
    
    img_dir = img_dirs{d};
    
    % collect images
    imgs = dir([img_dir '/*.nii']);
    imgs = {imgs.name};
    
    % get sub ids
    subs = {};
    for i = 1 : length(imgs)
        split_ = strsplit(imgs{i}, '_');
        subs{end+1} = split_{1};
    end
    subs = unique(subs);
    
    % average the images for each unique sub
    for s = 1 : length(subs)
        
        batch{1}.spm.util.imcalc.output = [subs{s} '_average+'];
        batch{1}.spm.util.imcalc.outdir = {img_dir};
        batch{1}.spm.util.imcalc.input = {[img_dir '/' subs{s} '_angle+.nii,1'];...
                                          [img_dir '/' subs{s} '_distance+.nii,1']};
        batch{1}.spm.util.imcalc.expression = '(i1+i2)/2';
        batch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        batch{1}.spm.util.imcalc.options.dmtx = 0;
        batch{1}.spm.util.imcalc.options.mask = 0;
        batch{1}.spm.util.imcalc.options.interp = 1;
        batch{1}.spm.util.imcalc.options.dtype = 4;
        
        spm_jobman('run', batch);
    end
end

