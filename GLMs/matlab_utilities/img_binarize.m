function img_binarize(img_path, thresh)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Binarize an image based on a threshold using SPM
    %
    % Arguments
    % ---------
    % img_path : str
    % thresh : numeric 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % minerva path, no harm no foul
    addpath /hpc/packages/minerva-centos7/spm/spm12
    
    [output_dir, img_fname, ~] = fileparts(img_path);
    splitname = strsplit(img_fname, '.nii');
    output_fname = [splitname{1} '_binarized-thresh' num2str(thresh * 100) '.nii'];
        
    spm('defaults', 'FMRI');
    batch{1}.spm.util.imcalc.input = cellstr(img_path);
    batch{1}.spm.util.imcalc.output = output_fname;
    batch{1}.spm.util.imcalc.outdir = {output_dir};
    batch{1}.spm.util.imcalc.expression = ['i1>' num2str(thresh)];
    batch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    batch{1}.spm.util.imcalc.options.dmtx = 0;
    batch{1}.spm.util.imcalc.options.mask = 0;
    batch{1}.spm.util.imcalc.options.interp = 1;
    batch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run', batch)