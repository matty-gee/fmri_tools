function cat_segmentation(t1_image)

    addpath /hpc/packages/minerva-centos7/spm/spm12

    batch{1}.spm.tools.cat.estwrite.data = {t1_image};
    batch{1}.spm.tools.cat.estwrite.data_wmh = {''};
    batch{1}.spm.tools.cat.estwrite.nproc = 2;
    batch{1}.spm.tools.cat.estwrite.useprior = '';
%     batch{1}.spm.tools.cat.estwrite.opts.tpm = {'/Users/matthew/Documents/MATLAB/spm12/tpm/TPM.nii'};
    batch{1}.spm.tools.cat.estwrite.opts.tpm = {'/hpc/packages/minerva-centos7/spm/spm12/tpm/TPM.nii'};
    batch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
    batch{1}.spm.tools.cat.estwrite.opts.biasacc = 0.5;

    batch{1}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.3];
    batch{1}.spm.tools.cat.estwrite.extopts.setCOM = 1;
    batch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
    batch{1}.spm.tools.cat.estwrite.extopts.affmod = 0;
    batch{1}.spm.tools.cat.estwrite.extopts.spm_kamap = 0;
    batch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
    batch{1}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
    batch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
    batch{1}.spm.tools.cat.estwrite.extopts.WMHC = 2;
%     batch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/Users/matthew/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
    batch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/hpc/packages/minerva-centos7/spm/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
    batch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
    batch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
    batch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
    batch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
    batch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
    batch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
    batch{1}.spm.tools.cat.estwrite.output.surface = 1;
    batch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;

    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
    batch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};

    batch{1}.spm.tools.cat.estwrite.output.GM.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
    batch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;

    batch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
    batch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;

    batch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
    batch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
    batch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
    batch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
    batch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.label.native = 1;
    batch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
    batch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
    batch{1}.spm.tools.cat.estwrite.output.las.native = 0;
    batch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
    batch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
    batch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
    batch{1}.spm.tools.cat.estwrite.output.warps = [1 0];
    batch{1}.spm.tools.cat.estwrite.output.rmat = 0;

    spm_jobman('run', batch);