function preprocess(sub_dir, sample)

if contains(pwd, 'sc/arion')
    addpath /sc/arion/projects/k23/code/matlab_utilities
    addpath /sc/arion/projects/k23/code/GLMs
    addpath /hpc/packages/minerva-centos7/spm/spm12
    spm_dir = '/hpc/packages/minerva-centos7/spm/spm12';
end

spm_jobman('initcfg');
spm('defaults', 'FMRI');

%% get data

func_imgs = cellstr(spm_select('ExtFPList', [sub_dir '/func/'], 'func.nii'));
anat_img  = cellstr(spm_select('FPList', [sub_dir '/anat/'], 'anat.nii')); 

%% set can acquisition parameters
% -- to read dicom use: dicominfo('dicom.dcm')

if sample == 1

    nslices     = 36;
    slice_order = [2:2:nslices 1:2:nslices]; %siemens even/odd ascending interleaved
    ref_slice   = slice_order(nslices/2); % middle of TR in slice ix
    vox         = [3 3 3];
    tr          = 2;
    ta          = tr-(tr/nslices);

elseif sample == 2

    nslices     = 70;
    slice_order = [0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5; % slice acquisition order is in ms for multiband data
                   0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5;
                   0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5;
                   0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5;
                   0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5;
                   0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5;
                   0; 685; 392.5; 97.5; 782.5; 490.00000001; 195.00000001; 882.50000001; 587.49999999; 292.5];
    ref_slice   = slice_order(nslices/2); % middle of TR in ms
    vox         = [2.1 2.1 2.1];
    tr          = 1;
    ta          = tr-(tr/nslices);

end


%% preprocessing

% STEPS:
% - 0: specify origin as anterior commissure 
% - 1: re-align & unwarp
% - 2: slice-time correction
% - 3: co-registration
% - 4: segmentation
% - 5: normalization
    
%--------------------------------------------------------------------------------------------------------------
% re-align & unwarp

% realignment: realign timeseries to ref image with least squares and 6 parameter
% - (rigid body) spatial transformation

% unwarping: predicated on idea that residual movement related variance
% - (variance explained by motion parameters, after realignment) is to
% - due to susceptibility by movement interaction creating deformations
% - in field 
%--------------------------------------------------------------------------------------------------------------

batch{1}.spm.spatial.realignunwarp.data.scans = func_imgs;
batch{1}.spm.spatial.realignunwarp.data.pmscan = '';
batch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
batch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
batch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
batch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
batch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
batch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
batch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
batch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
batch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
batch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
batch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
batch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
batch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
batch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
batch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
batch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
batch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
batch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
batch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
batch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
batch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
batch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

%--------------------------------------------------------------------------------------------------------------
% slice-time correction
%--------------------------------------------------------------------------------------------------------------

batch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
batch{2}.spm.temporal.st.nslices = nslices;
batch{2}.spm.temporal.st.tr = tr;
batch{2}.spm.temporal.st.ta = ta;
batch{2}.spm.temporal.st.so = slice_order;
batch{2}.spm.temporal.st.refslice = ref_slice;
batch{2}.spm.temporal.st.prefix = 'a';

%--------------------------------------------------------------------------------------------------------------
% co-registration: align t1 & functional images
% - some diffs from realignment: images prob have slightly diff shapes
% -- and diff tissue intensities & therefore have some mre steps

% - ref: reference image that remains stationary (mean image from realignment & unwarping)
% - source: image moved to match the ref (t1 image)
% - other: remain in alignment with source img (get same transformation as source) (realigned & slice time corrected t2*)
%--------------------------------------------------------------------------------------------------------------

batch{3}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
batch{3}.spm.spatial.coreg.estimate.source = anat_img;
batch{3}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
batch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi'; % normalized mutual information
batch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
batch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
batch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

%--------------------------------------------------------------------------------------------------------------
% segmentation: part of the unified segmentation & normalization algorithm in spm 

% tissue types: 1 = grey matter, 2 = white matter, 3 = csf

% options to save images:
% - batch{b}.spm.spatial.preproc.tissue(i).native: output tissue class image in native space: [native(c*) dartel_imported(rc*)]
% - batch{b}.spm.spatial.preproc.tissue(i).warped: output tissue class image in normalized space: [modulated(mwc*) unmodulated(wc*)]
% -- files can be used for vbm: smooth & estimate stats

% prefixes:
% - w = warped/normalized space
% - c = native/subject space
% - r = for dartel toolbox
% - m = modulated
% -- preserve gm signal in normalized partitions by compensating for volumetric changes induced by spatial normalization
%--------------------------------------------------------------------------------------------------------------

native = [1 1];
warped = [1 1];

batch{4}.spm.spatial.preproc.channel.vols(1) = anat_img;
batch{4}.spm.spatial.preproc.channel.biasreg = 0.001; % bias regularization: if data has v. little intensity non-uniformity artifact, use more severe regularization (alg. will try less hard to model it)
batch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
batch{4}.spm.spatial.preproc.channel.write = [0 1]; % save [bias_field bias_corrected]

batch{4}.spm.spatial.preproc.tissue(1).tpm = {[spm_dir '/tpm/TPM.nii,1']};
batch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
batch{4}.spm.spatial.preproc.tissue(1).native = native; 
batch{4}.spm.spatial.preproc.tissue(1).warped = warped; 

batch{4}.spm.spatial.preproc.tissue(2).tpm = {[spm_dir '/tpm/TPM.nii,2']};
batch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
batch{4}.spm.spatial.preproc.tissue(2).native = native;
batch{4}.spm.spatial.preproc.tissue(2).warped = warped;

batch{4}.spm.spatial.preproc.tissue(3).tpm = {[spm_dir '/tpm/TPM.nii,3']};
batch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
batch{4}.spm.spatial.preproc.tissue(3).native = native;
batch{4}.spm.spatial.preproc.tissue(3).warped = warped;

batch{4}.spm.spatial.preproc.tissue(4).tpm = {[spm_dir '/tpm/TPM.nii,4']};
batch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
batch{4}.spm.spatial.preproc.tissue(4).native = native;
batch{4}.spm.spatial.preproc.tissue(4).warped = warped;

batch{4}.spm.spatial.preproc.tissue(5).tpm = {[spm_dir '/tpm/TPM.nii,5']};
batch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
batch{4}.spm.spatial.preproc.tissue(5).native = native;
batch{4}.spm.spatial.preproc.tissue(5).warped = warped;

batch{4}.spm.spatial.preproc.tissue(6).tpm = {[spm_dir '/tpm/TPM.nii,6']};
batch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
batch{4}.spm.spatial.preproc.tissue(6).native = native;
batch{4}.spm.spatial.preproc.tissue(6).warped = warped;

% get segmentation forward deformation field that allows normalization to mni space
batch{4}.spm.spatial.preproc.warp.mrf = 1;
batch{4}.spm.spatial.preproc.warp.cleanup = 1; % clean up: routine for extracting brain from segmented images; 1 = light clean up, 2 = thorough clean up
batch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]; % warping regularization
batch{4}.spm.spatial.preproc.warp.affreg = 'mni'; % affine regularization - IMPORTANT: images should be first be placed in approx. alignment with MNI using display function in SPM
batch{4}.spm.spatial.preproc.warp.fwhm = 0; % "usually can have 0 for fmri"
batch{4}.spm.spatial.preproc.warp.samp = 3;
batch{4}.spm.spatial.preproc.warp.write = [1 1]; % save forward (to normalize to MNI) & inverse (to normalize, eg, GIFTI files) deformation fields: [inverse forward]

%--------------------------------------------------------------------------------------------------------------
% normalization: apply segmentation forward deformation on coregistered imgs
%--------------------------------------------------------------------------------------------------------------

batch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
batch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
batch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
batch{5}.spm.spatial.normalise.write.woptions.vox = vox; % can just be voxel size: can't achieve higher resolution than native dimensions, so smaller values will just increase size of files
batch{5}.spm.spatial.normalise.write.woptions.interp = 4; % interpolation default = 4th degree B-Spline; higher degrees are slower b/c use more neighbors
batch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';

%--------------------------------------------------------------------------------------------------------------
% run it
%--------------------------------------------------------------------------------------------------------------

spm_jobman('run', batch);

