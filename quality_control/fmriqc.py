#!/usr/bin/python

"""
    modified from pyfmri qc
    options
        -n:   functional MR nifti file
        -s:   percentage of voxel with the lowest values (outside the mask) for SNR calculation
        either
            -t:   threshold of minimum mean values of voxel that should be included in the quality check
        or
            -k:   binary nifti mask file of voxel that should be included in the quality check
        optional
        -o:   output directory
        -m:   motion parameters file of motion correction from FSL (*.par), SPM (rp*.txt) or AFNI (*.1D).

    Example:
        python pyfMRIqc.py

    OUTPUT
    (all output files ends with the input filename before the file extension)

    png images:
        pyfMRIqc_ ... .py:
        - MEAN_<yourfile>
            showing the mean voxel intensity over time in axial slices
        - VARIANCE(_<yourfile>_thr<XXX>)
            showing variance of the voxel time courses (max threshold XXX) in axial slices
        - MASK_<yourfile>
            showing the voxels included in the QC (based on -k mask or -t threshold input) in blue and the n% of voxels with
            lowest mean voxel intensity (based on -s input) in green on top of the mean image.
        - BINMEAN_<yourfile>
            The value range of mean voxel intenstiy is devided in 50 bins with equal number of voxels.
            This image shows the average time course over voxels for each of the bins.
        - SUM_SQUARED_SCALED_DIFF_<yourfile>
        - PLOTS_<yourfile> containing the following:
            - scaled variability: Mean (over all voxel) squared diff plot over time / global mean squared diff
            - slice by slice variability: Mean (mean per slice) squared diff plot over time / global mean squared diff
            - if motion file was selected: sum of relative movements over time (z-scored)
            - scaled mean voxel intensity: mean(data/global mean) )(z-scored)
            - variance of scaled variability (z-scored)
            - min/mean/max slice variability

    html file
        pyfMRIqc_HTML_<yourfile>.html containing:
        - Overview about scan and QC parameters
        - All png images explained above
        - Summary of signal to noise ratio (SNR) calculation
        - Summary of motion parameter (if motion parameter file was specified as input)

    text file
        pyfMRIqc_textfile_<yourfile>.txt containing an overview about scan, QC and motion parameters (similar to the html
        file)

    nifti images:
        - mean_<yourfile>
            mean voxel intensity over time (3D)
        - variance_<yourfile>
            variance of the voxel time courses(3D)
        - mask_<yourfile>
            binary - containing voxels above the threshold or the input mask (3D)
        - mask4snr_<yourfile>
            binary - lowest n percent of lowest values used for SNR calculation (3D)
        - snr_<yourfile>
            voxel-wise signal-to-noise ratio (3D)
        - squared_scale_diff_<yourfile>
            squared scaled signal variability: squared difference between two consecutive volumes divided by the global
            mean difference
"""
# matthew schafer

import os, sys, copy, glob
import nibabel as nib
import numpy as np
import pandas as pd
from datetime import datetime
from scipy import stats
# import easygui
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
from nilearn import plotting
from nilearn.maskers import NiftiSpheresMasker, NiftiMasker

#----------------------------------------------------------
# motion specific
#----------------------------------------------------------

def calc_fd(movement_regressors, volumes_back=1):

    ''' Calculate framewise displacement with respect to a specified number of volumes back 
        fast TR fMRI might need to account for diff. kinds of motions etc...
        https://www.sciencedirect.com/science/article/pii/S1053811919306226?via%3Dihub#fig5
    '''

    assert movement_regressors.shape[1] == 6, f'movement_regressors shape is {movement_regressors.shape}'
    if isinstance(movement_regressors, pd.DataFrame):
        movement_regressors = movement_regressors.values
    fd = np.zeros(movement_regressors.shape[0])

    # Calculate FD, starting from the index corresponding to volumes_back
    for i in range(volumes_back, movement_regressors.shape[0]):
        # Calculate the sum of absolute differences between the current and the volumes_back-th previous volume
        fd[i] = np.sum(np.abs(movement_regressors[i,:] - movement_regressors[i - volumes_back,:]))

    fd[:volumes_back] = np.nan # Handling the initial volumes by assigning np.nan
    return fd[:, np.newaxis]

# def calc_fd(rp):
#     ''' calculate framewise displacement from motion parameters '''
#     # expects rp to have shape (num_trs, 6)
#     assert rp.shape[1] == 6, f'rp shape is {rp.shape}'
#     fd = np.sum(np.abs(np.diff(rp, axis=0)), axis=1) # sum of absolute value displacement from previous volume
#     fd = np.hstack([np.nan, fd]) # bc 1st vol has no previous vol to compare to
#     assert fd.shape == (rp.shape[0],), f'fd shape is {fd.shape}'
#     return fd[:,np.newaxis]

def save_rp_plot(rp_fname, windows=None, out_dir=None):

    # windows: list of lists, each list is a window of TRs to highlight
    # for social space project specifically
    sub_id = rp_fname.split('/')[-3]

    # load in rps, add fd
    try: 
        rp = pd.read_csv(rp_fname, header=None, delim_whitespace=True)
        rp.columns = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']
    except: 
        rp = pd.read_csv(rp_fname)
        rp.columns = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']
    rp.insert(0, 'fd', calc_fd(rp))

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(30, 8))
    sns.lineplot(data=rp)
    if windows is not None:
        for window in windows:
            ax.axvspan(window[0], window[-1], color='grey', alpha=0.3, label='decision')
    ax.set_title('Realignment Parameters', fontsize=20)
    ax.set_xlabel('TR', fontsize=16)
    ax.set_ylabel('mm/degrees', fontsize=16)

    min_ = np.min(rp[['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']].values) - 0.5
    max_ = np.max([np.max(rp['fd']), 2.5]) # if 2.5 is bigger then use that
    ax.set_ylim(-max_, max_)

    # add a hortixontal line for each 1 mm
    for i in range(1, int(max_)+1):
        ax.axhline(i, color='black', linestyle='--', alpha=0.5)
        ax.axhline(-i, color='black', linestyle='--', alpha=0.5)

    # # annotate with the max and mean of each in rp
    # for i, col in enumerate(rp.columns):
    #     max_ = np.max(rp[col])
    #     mean_ = np.mean(rp[col])
    #     ax.annotate(f'{col} max={max_:.2f}, mean={max_:.2f}', 
    #                 xy=(0.01, 0.98-(i*.05)), xycoords='axes fraction', fontsize=13)

    # only show the first 8 legend items
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:8], labels[:8], loc='center left', bbox_to_anchor=(0.0, 0.18), frameon=True)

    # savefig
    if out_dir is None:
        out_dir = f'{os.path.dirname(rp_fname)}/QC' 
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    fig.savefig(f'{out_dir}/{sub_id}_motion_plot.png', bbox_inches='tight')
    plt.close()


#----------------------------------------------------------
# SNR images
#----------------------------------------------------------


def process_qc(nii_fname, rp_fname, mask_nii_fname=None, 
            mask_thresh=None, snr_voxel_percentage=None, 
            out_dir=None, plot=False):

    print('Running QC for fMRI data')

    #--------------------------------------------------
    # clean up inputs
    #--------------------------------------------------
    
    filepath, filename = os.path.split(nii_fname)
    x, fx = os.path.splitext(filename)
    if fx == ".gz":
        fname, _ = os.path.splitext(x)
        fext = ".nii.gz"
    else:
        fname = x
        fext = ".nii"

    if not out_dir: 
        out_dir = os.path.join(filepath, "QC")
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    print(f'Functional MRI nifti file: {nii_fname}')
    print(f'Motion parameter file: {rp_fname}')
    if mask_thresh is not None: 
        print(f'Mask threshold: {str(mask_thresh)}')
    if mask_nii_fname is not None: 
        print(f'Mask nifti file: {mask_nii_fname}')
    if snr_voxel_percentage is not None:
        print(f'SSNR voxel percentage: {str(snr_voxel_percentage)}')
    else:
        print('SSNR based on non-brain ROI')
    print(f'Output directory: {out_dir}')

    #--------------------------------------------------
    # Load and get func data
    #--------------------------------------------------

    print("Loading functional image")
    nii    = nib.load(nii_fname)
    header = nii.header
    affine = nii.affine
    voxelsize = header['pixdim'][1]
    data      = nii.get_fdata()
    shape     = np.array(data)[:, :, :, 0].shape
    nrvoxel   = shape[0] * shape[1] * shape[2]
    nrvolumes = np.array(data)[0, 0, 0, :].shape
    del nii

    #--------------------------------------------------
    # Create mean & std data & nifti images
    #--------------------------------------------------

    meandata = np.mean(data, axis=3) # mean over time
    print(f'Creating mean image, w/ shape={meandata.shape}')
    header2 = header
    header2['glmin'], header2['glmax'] = np.min(meandata), np.max(meandata)
    mean_img = nib.Nifti1Image(meandata, affine, header2)
    nib.save(mean_img, os.path.join(out_dir, f'mean_{fname}{fext}'))

    stddata = np.std(data, axis=3) # mean over time
    print(f'Creating std image, w/ shape={stddata.shape}')
    header2 = header
    header2['glmin'], header2['glmax'] = np.min(stddata), np.max(stddata)
    std_img = nib.Nifti1Image(stddata, affine, header2)
    nib.save(std_img, os.path.join(out_dir, f'std_{fname}{fext}'))

    # create png
    if plot:
        meanimage = nii2image(meandata, 'Mean', os.path.join(out_dir, f'MEAN_{fname}.png'))
        meanimage = meanimage / np.max(meanimage[:]) * 255
        nrbins = 50

        # bins equal number of voxel
        binmean = np.zeros((nrbins, nrvolumes[0]))
        d = meandata.reshape((nrvoxel, ))
        di = np.argsort(d)
        dbins = np.linspace(0, nrvoxel, nrbins)
        dbinlabs = np.zeros((nrbins,))

        for nn in range(nrbins):
            m = np.zeros((nrvoxel, 1))
            try:
                m[int(np.round(dbins[nn], decimals=0)):int(np.round(dbins[nn+1], decimals=0))] = 1
            except:
                m[int(np.round(dbins[nn], decimals=0)):] = 1
            binmask = m[di.argsort()]
            binmask = binmask.reshape((shape[0], shape[1], shape[2]))
            k = np.where(binmask == 1)
            if len(k[0]) > 0:
                binmean[nn, :] = stats.zscore(np.mean(data[k[0][:], k[1][:], k[2][:], :], axis=0))
                dbinlabs[nn] = np.max(meandata[k[0][:], k[1][:], k[2][:]])
            del k
            del binmask
            del m

        binticks = np.arange(1, nrbins, 5)
        binticklabels = dbinlabs[1::5]
        binticklabels = binticklabels.astype(int)

        plt.figure(num=None, figsize=(5, 3.5), dpi=300, facecolor=(210 / 255, 227 / 255, 244 / 255), edgecolor='k')
        plt.imshow(binmean, cmap='gray', aspect='auto')
        plt.xlabel("volumes")
        plt.ylabel("bins")
        plt.yticks(binticks, binticklabels, fontsize=5)
        plt.xticks(fontsize=5)
        plt.title("bins with equal number of voxel (" + str(int(nrvoxel/nrbins)) + " per bin)")
        plt.savefig(os.path.join(out_dir, f'BINMEAN_{fname}.png'), dpi=300, facecolor=(210 / 255, 227 / 255, 244 / 255))

    #--------------------------------------------------
    # Create SNR images
    #--------------------------------------------------
    # TODO add a temporal contrast snr measure, with the decision trials vs. the non-decision trials
    # http://www.newbi4fmri.com/mini-tutorial-signal-to-noise


    # # noise mask: voxels w/ smallest signal ... but here some of the voxels are negative?
    # snrvoxel_n = int(nrvoxel * snr_voxel_percentage / 100)
    # vec = np.sort(meandata.reshape(1, nrvoxel), axis=1) # flatten & sort
    # snrvoxel_mean, snrvoxel_std = np.mean(vec[0][0:snrvoxel_n]), np.std(vec[0][0:snrvoxel_n])
    # snrvoxel_min, snrvoxel_max  = vec[0][0], vec[0][snrvoxel_n]
    # ssnrmask = np.where((meandata > 0) & (meandata < snrvoxel_max), 1, 0) # voxels w/in noise range
    # stdnoise = np.std(meandata[ssnrmask == 1]) # std over noisiest voxels
    # print(f'Noise min={snrvoxel_min}, max={snrvoxel_max}, mean={snrvoxel_mean}, std={snrvoxel_std}')
    # print(f'{np.sum(ssnrmask)} voxels included in noise mask')
    # nib.save(nib.Nifti1Image(ssnrmask, affine, header), 
    #          os.path.join(out_dir, f"mask4ssnr_{fname}{fext}"))
    
    print("Calculating spatial SNR image")
    if snr_voxel_percentage is None: # non-brain ROI (5mm outside brain)
        try: 
            adj = np.round(18 / voxelsize) # 18mm from top right corner - should have enough space for 5mm radius  
            nonbrain_coords = shape - np.array((adj, adj, adj)) 
            ssnr_masker = NiftiSpheresMasker([nonbrain_coords], radius=5)
            noise = ssnr_masker.fit_transform(nii_fname)
        except:
            adj = np.round(25 / voxelsize) # 18mm from top right corner - should have enough space for 5mm radius  
            nonbrain_coords = shape - np.array((adj, adj, adj))
            ssnr_masker = NiftiSpheresMasker([nonbrain_coords], radius=5)
            noise = ssnr_masker.fit_transform(nii_fname)
        print(f'Non-brain ROI xyz: {nonbrain_coords}')
        snrvoxel_mean, snrvoxel_std = np.mean(noise), np.std(noise)
        snrvoxel_min, snrvoxel_max = np.min(noise), np.max(noise)
        snrvoxel_n = len(noise)
        stdnoise = np.std(noise)
        print(f'{snrvoxel_n} voxels included in noise mask')
    ssnr_data = np.divide(meandata, stdnoise).astype(np.int16)
    header2 = header
    header2['glmin'], header2['glmax'] = np.min(ssnr_data), np.max(ssnr_data)
    nib.save(nib.Nifti1Image(ssnr_data, affine, header2), 
             os.path.join(out_dir, f'ssnr_{fname}{fext}'))
    
    print("Calculating temporal SNR image") # timeseries mean / timeseries std
    tssnr_data = np.divide(meandata, stddata).astype(np.int16)
    header2  = header
    header2['glmin'], header2['glmax'] = np.min(tssnr_data), np.max(tssnr_data)
    nib.save(nib.Nifti1Image(tssnr_data, affine, header2), 
            os.path.join(out_dir, f'tsnr_{fname}{fext}'))

    # plot the images
    for snr_type in ['tsnr', 'ssnr']:
        plotting.plot_stat_map(f'{snr_type}_{fname}{fext}', display_mode='mosaic',
                            title="display_mode='mosaic' default cut_coords")
        plt.savefig(f'{out_dir}/{snr_type}_{fname}_slices.png')
      
    # brain mask: voxels to include in quality check
    mask = []
    if mask_thresh is not None:

        # check overlap between mask_thresh and snr max value
        if snrvoxel_max > mask_thresh:
            print("INPUT ERROR: Given mask threshold is smaller than max intensity value of given percentage of "
                  "voxel for SNR calculation.\n"
                  "SNR value range of " + str(snr_voxel_percentage) +
                  "% voxel with lowest intensity = 0 - " + str(snrvoxel_max) + "\n"
                  "Total value range = 0 - " + str(np.max(data[:])) + "\n"
                  "Specified mask threshold = " + str(mask_thresh) + "\n"
                  "Either increase threshold or decrease percentage of voxel for SNR.")
            sys.exit(2)
        
        print("Create MASK")
        mask = np.where(meandata >= mask_thresh, 1, 0)
    
    elif mask_nii_fname is not None:
        
        print("Load Nifti Mask")
        masknii = nib.load(mask_nii_fname)
        maskdata = masknii.get_fdata()
        mask = np.asarray(maskdata)
        if len(maskdata.shape) == 4:
            maskshape = np.array(maskdata)[:, :, :, 0].shape
            mask = mask[:, :, :, 0]
        elif len(maskdata.shape) == 3:
            maskshape = np.array(maskdata)[:, :, :].shape
            mask = mask[:, :, :]
        del masknii

        print("Check Nifti Mask")
        maskunique = np.unique(maskdata)
        if np.min(maskunique) < 0:
            easygui.msgbox("Warning", "Mask contains value < 0")
        elif np.max(maskunique) > 1:
            easygui.msgbox("Warning", "Mask contains value > 1")
        elif len(maskunique) != 2:
            easygui.msgbox("Warning", "Mask is not binary")
        elif sum(np.asarray(shape) - np.asarray(maskshape)) != 0:
            easygui.msgbox("Mask dimensions not equal to functional data")

        if mask_thresh is not None:
            if np.min(meandata[mask == 1]) > mask_thresh:
                print("INPUT ERROR: Mask contains voxel with smaller intensity than max intensity value of given "
                    "percentage of voxel for SNR calculation.\n"
                    "SNR value range of " + str(snr_voxel_percentage) +
                    "% voxel with lowest intensity = 0 - " + str(snrvoxel_max) + "\n"
                    "Total value range = 0 - " + str(np.max(data[:])) + "\n"
                    "Minimum voxel intensity in mask = " + str(np.min(meandata[mask == 1])) + "\n"
                    "Either change mask or decrease percentage of voxel for SNR.")
                sys.exit(2)
    
    nib.save(nib.Nifti1Image(mask, affine, header), 
             os.path.join(out_dir, f"mask_{fname}{fext}"))

    # apply mask for downstream operations
    s = np.asarray(np.asarray(np.array(data)).shape)
    for ii in range(s[3]):
        data[:, :, :, ii] = np.multiply(mask, data[:, :, :, ii])

    # create pngs restricted to quality check region   
    if plot: 
        pngfilename = os.path.join(out_dir, f'MASK_{fname}.png')
        maskimage = nii2image(mask, 'Mask', pngfilename)
        imageshape = np.array(meanimage).shape
        meanimagecol = np.zeros((imageshape[0], imageshape[1], 3))
        meanimagecol[:, :, 0] = meanimage

        # add SNR mask in green
        pngfilename = os.path.join(out_dir, f'SNR_{fname}.png')
        snrimage = nii2image(ssnrmask, 'SNR', pngfilename)
        meanimage0 = copy.deepcopy(meanimage)
        meanimage0[snrimage > 0] = 255
        meanimagecol[:, :, 1] = meanimage0

        # add mask in blue
        meanimage2 = copy.deepcopy(meanimage)
        meanimage2[maskimage > 0] = 255
        meanimagecol[:, :, 2] = meanimage2
        pngfilename = os.path.join(out_dir, f'MASK_{fname}.png')
        sizes = np.shape(meanimagecol)
        height = float(sizes[0])
        width = float(sizes[1])
        
        # plot
        fig = plt.figure()
        fig.set_size_inches(width / height, 1, forward=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(meanimagecol / 255.0)
        plt.savefig(pngfilename, dpi=height)

    #--------------------------------------------------
    # text file w/ summary info
    #--------------------------------------------------

    # open text file
    print("Creating text file")
    textfilename = os.path.join(out_dir, f"qc-summary_{fname}.txt")
    text_file = open(textfilename, "w")
    text_file.write("SCAN PARAMETERS AND USER INPUT: \n\n")
    text_file.write(f"Slice Resolution: {header['dim'][1]}x{header['dim'][2]}\n")
    text_file.write(f"Num slices: {header['dim'][3]}\n")
    text_file.write(f"Num Volumes: {data.shape[3]}\n")
    text_file.write(f"Voxel Size: {header['pixdim'][1]}x{header['pixdim'][2]}x{header['pixdim'][3]}\n")
    text_file.write(f"Total Num Voxels: {data.shape[0] * data.shape[1] * data.shape[2]}\n")
    text_file.write(f"Mask Threshold: {mask_thresh}\n")
    text_file.write(f"Num Mask Voxels: {np.sum(mask)}\n")
    if snr_voxel_percentage is not None: # use a voxel threshold
        text_file.write(f"Percentage of voxel with lowest values for sSNR: {snr_voxel_percentage}\n")
        text_file.write(f" sSNR nr of lowest voxel: {snrvoxel_n}\n")
    else: # created an ROI instead
        text_file.write(f"Created non-brain ROI for sSNR: {nonbrain_coords}\n")
    text_file.write(f" sSNR voxel MEAN: {snrvoxel_mean}\n")
    text_file.write(f" sSNR voxel STD: {snrvoxel_std}\n")
    text_file.write(f" sSNR voxel value range: {snrvoxel_min} - {snrvoxel_max}\n")
    if mask_nii_fname is not None:
        text_file.write(f"Mask file: {mask_nii_fname}\n")
    text_file.write("\n------------------------------------ \n\n")
    text_file.write("QC OVERVIEW: \n\n")
    text_file.write(f"Mean: {np.mean(data)}\n")
    text_file.write(f"Mean (mask): {np.mean(data[mask == 1])}\n")
    text_file.write(f"SD: {np.std(data)}\n")
    text_file.write(f"SD (mask): {np.std(data[mask == 1])}\n")

    # add mean SNR over slice to text file
    snrvec = np.zeros((ssnr_data.shape[2], 1))
    for ii in range(ssnr_data.shape[2]):
        snrslice = ssnr_data[:, :, ii]
        maskslice = mask[:, :, ii]
        if np.sum(maskslice) > 0:
            snrvec[ii] = np.nanmean(snrslice[maskslice == 1])
        else:
            snrvec[ii] = np.nan
    text_file.write(f"\nSignal-To-Noise Ratio: \n")
    text_file.write(f"Min Slice SNR: {np.nanmin(snrvec)}\n")
    text_file.write(f"Max Slice SNR: {np.nanmax(snrvec)}\n")
    text_file.write("ALL Slice SNRs: ")
    for ii in range(len(snrvec)):
        text_file.write(f"{snrvec[ii]} ")
    text_file.write("\n")
    text_file.write(f"Mean voxel SNR: {np.nanmean(snrvec)}\n")

    #--------------------------------------------------
    # Process motion data (if specified)
    #--------------------------------------------------

    if rp_fname is not None:

        print("Loading motion parameters")

        relrm, rmsum, nrrm01, nrrm05, nrrmv = [], [], [], [], []
        x, motionext = os.path.splitext(rp_fname)
        rm = pd.read_csv(rp_fname, header=None, delim_whitespace=True)

        # motion in mm assuming head radius is 5cm
        if motionext == ".txt":  # SPM
            for ii in [3, 4, 5]:
                rm[:, ii] = rm[:, ii] * 50

        elif motionext == ".par":  # FSL
            for ii in [0, 1, 2]:
                rm[:, ii] = rm[:, ii] * 50

        elif motionext == ".1D":  # AFNI
            for ii in [0, 1, 2]:
                rm[:, ii] = np.radians(rm[:,ii]) # convert from degree to radians
                rm[:, ii] = rm[:, ii] * 50

        relrm = np.diff(rm, axis=0) # create relative values
        relrm = np.absolute(relrm) # create absolute values of relative movement
        rm = np.absolute(rm) # absolute values of absolute movement
        rmsum = np.sum(rm, axis=1) # sum absolute values

        # get thresholds:
        nrrm01 = relrm[np.where(relrm > .1)]
        nrrm05 = relrm[np.where(relrm > .5)]
        nrrmv  = relrm[np.where(relrm > voxelsize)]
   
        # add movent to textfile
        text_file.write("\nMovement: \n")
        text_file.write("Mean absolute Movement: " + str(np.mean(rm)) + "\n")
        text_file.write("Max absolute Movement: " + str(np.max(rm)) + "\n")
        text_file.write("Mean relative Movement: " + str(np.mean(relrm)) + "\n")
        text_file.write("Max relative Movement: " + str(np.max(relrm)) + "\n")
        text_file.write("Relative movements (>0.1mm): " + str(len(nrrm01)) + "\n")
        text_file.write("Relative movements (>0.5mm): " + str(len(nrrm05)) + "\n")
        text_file.write("Relative movements (>voxelsize): " + str(len(nrrmv)) + "\n")

    text_file.close()

    #----------------------------------------------------- 
    # variance
    #-----------------------------------------------------

    print("Calculating variance")
    vardata = np.var(data, axis=3).astype(np.int32)
    header2 = header
    header2['glmin'], header2['glmax'] = np.nanmin(vardata), np.nanmax(vardata)
    nib.save(nib.Nifti1Image(vardata, affine, header2), 
             os.path.join(out_dir, f"variance_{fname}{fext}"))
    
    del vardata

    if plot:
        varfilename = os.path.join(out_dir, 'VARIANCE_' + fname + '.png')
        varthresh = nii2image(vardata, 'Variance', varfilename)
    
    #-----------------------------------------------------
    # scaled squared diff
    #-----------------------------------------------------

    print("Calculating scaled squared difference")
    diff2data = np.power(np.diff(data, axis=3), 2).astype(np.int32)
    diff2data_a = np.divide(diff2data, np.mean(diff2data)).astype(np.int32)
    del diff2data
    header2 = header
    header2['glmin'], header2['glmax'] = np.min(diff2data_a), np.max(diff2data_a)
    nib.save(nib.Nifti1Image(diff2data_a, affine, header2), 
            os.path.join(out_dir, "squared_diff_scaled_" + fname + fext))

    # create png
    if plot:
        sumdiff2data_a = np.sum(diff2data_a, axis=3)
        pngfilename = os.path.join(out_dir, 'SUM_SQUARED_DIFF_SCALED_' + fname + '.png')
        mssdthresh = nii2image(sumdiff2data_a, 'DIFF', pngfilename)

        # plot data
        print("Creating plot")
        plotnr = 4
        fig = plt.figure(frameon=False)
        fig.set_size_inches(12, 14)

        # plot 1
        d1 = np.mean(np.mean(np.mean(diff2data_a, axis=0), axis=0), axis=0)
        plt.subplot(plotnr, 1, 1)
        plt.plot(d1)
        plt.xlabel('Difference image number')
        plt.ylabel('mean SSD')

        # plot 2
        d2 = np.mean(np.mean(diff2data_a, axis=0), axis=0)
        plt.subplot(plotnr, 1, 2)
        plt.plot(d2.T, 'x')
        plt.xlabel('Difference image number')
        plt.ylabel('slice-wise mean SSD')

        # plot 3
        d3 = np.mean(np.mean(np.mean(np.divide(data, np.mean(data)), axis=0), axis=0), axis=0)
        vardiff = np.zeros(diff2data_a.shape[3])
        for ii in range(diff2data_a.shape[3]):
            vardiff[ii] = np.var(diff2data_a[:, :, :, ii])
        plt.subplot(plotnr, 1, 3)
        plt.plot(stats.zscore(vardiff), label='variance of SSD')
        plt.plot(stats.zscore(d3), label='normalized mean voxel intensity')
        plt.xlabel('Image number')
        plt.ylabel('Normalized Amplitudes')
    
        if rp_fname is not None:

            # add line to plot
            relrmsum = np.sum(relrm, axis=1)
            plt.plot(stats.zscore(relrmsum), label='sum of relative movements')
            plt.plot(stats.zscore(rmsum), label='sum of absolute movements')

        plt.legend(loc='lower right')

        # plot 4
        d4a = np.max(d2, axis=1)
        d4b = np.min(d2, axis=1)
        d4c = np.mean(d2, axis=1)
        plt.subplot(plotnr, 1, 4)
        plt.plot(d4b, label='min')
        plt.plot(d4a, label='max')
        plt.plot(d4c, label='mean')
        plt.xlabel('Slice number')
        plt.ylabel('min/mean/max slice SSD')
        plt.legend(loc='upper right')

        # save and show figure
        newfilename = os.path.join(out_dir, "PLOTS_" + fname + ".png")
        plt.savefig(newfilename, facecolor=(210 / 255, 227 / 255, 244 / 255))

    #-----------------------------------------------------
    # create html file
    #-----------------------------------------------------

    # Open html file
    print("Create html file")
    htmlfilename = os.path.join(out_dir, "HTML_" + fname + ".html")
    html_output = open(htmlfilename, 'w')

    # add head to html file
    html_output.write("<html><head><body style=""background-color:#d2e3f4;""><title>pyfMRIqc output</title>")
    bootstraplines = """<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js" 
        integrity="sha384-ChfqqxuZUCnJSK3+MXmPNIyE6ZbWh2IMqE241rYiqJxyMiZ6OW/JmZQ5stwEULTy" 
        crossorigin="anonymous"></script> <link rel="stylesheet" 
        href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" 
        integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">"""
    html_output.write(bootstraplines)
    
    # add style for tables to html file
    style_text = """
        <style>
        table, th, td {
            border: 1px solid black;
            border-collapse: collapse;
            width: 400;    
            background-color: #fff;
            align="center"
        }
        th, td {
            padding: 15px;
            text-align: left;
        }
        table th {
            background-color: #eee;
        }
        hr {
            display: block;
            height: 1px;
            border: 20;
            border-top: 1px 
            solid #333;
            margin: 1em 0;
            padding: 0; 
        }
        .center {
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 90%;
        }
        h1 { padding-top: 45px; }
        h2 { padding-top: 20px; }
        h3 { padding-top: 25px; }
        body{font-size:20px;}
        </style>
        </head>
        <body>
    """
    html_output.write(style_text)

    navbar = """
            <nav class="navbar fixed-top navbar-expand-lg navbar-light bg-light">
            <a class="navbar-brand" href="#">pyfMRIqc</a>
            <div class="collapse navbar-collapse">
                <ul class="navbar-nav">
                <li class="nav-item"><a class="nav-link" href="#Input parameter">Input parameter</a></li>
                <li class="nav-item"><a class="nav-link" href="#Scan parameter">Scan parameter</a></li>
                <li class="nav-item"><a class="nav-link" href="#QC plots">QC plots</a></li>
                <li class="nav-item"><a class="nav-link" href="#BINMEAN">BINMEAN</a></li>
                <li class="nav-item"><a class="nav-link" href="#Mean">Mean</a></li>
                <li class="nav-item"><a class="nav-link" href="#Masks">Masks</a></li>
                <li class="nav-item"><a class="nav-link" href="#Variance">Variance</a></li>
                <li class="nav-item"><a class="nav-link" href="#SNR">SNR</a></li>
                <li class="nav-item"><a class="nav-link" href="#SUMDIFF">SUMDIFF</a></li>
                <li class="nav-item"><a class="nav-link" href="#Motion">Motion</a></li>
                <li class="nav-item"><a class="nav-link" href="#About">About</a></li>
                </ul>
            </div>
            </nav>"""
    html_output.write(navbar)

    # add parameter table to html file
    html_output.write("""<div id="Input parameter"> <h1>Input parameter</h1>""")
    parameter_table = """
        <table>
            <tr>
                <th>Functional file</th>
                <td>""" + fname + """</td>
            </tr>
            <tr>
                <th>Motion file</th>
                <td>""" + str(rp_fname) + """</td>
            </tr>
            <tr>
                <th>Threshold value</th>
                <td>""" + str(mask_thresh) + """</td>
            </tr>
            <tr>
                <th>Mask file</th>
                <td>""" + str(mask_nii_fname) + """</td>
            </tr>
            <tr>
                <th>SNR threshold</th>
                <td>""" + str(snr_voxel_percentage) + """</td>
            </tr>
        </table>"""
    html_output.write(parameter_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add scan parameters and user input to html file
    html_output.write("""<div id="Scan parameter"> <h1>Scan parameter</h1>""")
    acquisition_table = """
        <table>
            <tr>
                <th>Slice Dimensions</th>
                <td>""" + str(header['dim'][1]) + "x" + str(header['dim'][2]) + """</td>
            </tr>
            <tr>
                <th>Number of Slices</th>
                <td>""" + str(header['dim'][3]) + """</td>
            </tr>
            <tr>
                <th>Number of Volumes</th>
                <td>""" + str(data.shape[3]) + """</td>
            </tr>
            <tr>
                <th>Voxel Size</th>
                <td>""" + str(header['pixdim'][1]) + "x" + str(header['pixdim'][2]) + "x" + str(header['pixdim'][3]) +\
                        """</td>
            </tr>
            <tr>
                <th>Total Number of Voxels</th>
                <td>""" + str(data.shape[0] * data.shape[1] * data.shape[2]) + """</td>
            </tr>
        
        </table>"""
    html_output.write(acquisition_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # add QC plot
    html_output.write("""<div id="QC plots"> <h1>QC plots</h1>""")

    plottext = """<p>The first plot shows the mean of the scaled squared difference (SSD)over all voxel in 
                    the mask. The second plot shows the mean of the (SSD) for each slice separately. The third plot 
                    shows: 1) Normalised average of the demeaned voxel intensity of each volume. 2) Normalised 
                    variance of the SSD over all voxels in the mask. 3) (only if motion file was specified as input)
                    : Normalized sum of relative movement (sum of all relative motion translations and rotations.).
                    See <a href = "https://drmichaellindner.github.io/pyfMRIqc/#qc-plots">here</a>
                    for guidelines about how to use these plots.</p>"""
    html_output.write(plottext)
    html_output.write("""<img src ="pyfMRIqc_PLOTS_""" + fname + """.png" alt="pyfMRIqc plots" class="center">""")

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add BIN x VOLUME to html file
    html_output.write("""<div id="BINMEAN"> <h1>Mean voxel time course of bins with equal number of voxels</h1>""")
    bintext = """<p>Mean voxel time course of bins.  
                        See <a href = "https://drmichaellindner.github.io/pyfMRIqc/#mean-voxel-time-course-of-bins-with-equal-number-of-voxels">here</a>
                        for guidelines about how to use this plot.</p>"""
    html_output.write(bintext)
    html_output.write("""<img src="pyfMRIqc_BINMEAN_""" + fname +
                      """.png" alt="Mean signal from functional image" class="center">""")
    html_output.write("</div><br><hr><br>")

    # Add mean data to html file
    html_output.write("""<div id="Mean"> <h1>Mean voxel intensity</h1>""")
    html_output.write("""<img src="pyfMRIqc_MEAN_""" + fname +
                      """.png" alt="Mean signal from functional image" class="center">""")
    html_output.write("<h2> Mean voxel intensity summary</h2>")
    mean_signal_table = """
        <table style="width:50%">
            <tr>
                <th>Mean Signal (unmasked)</th>
                <td>""" + str(np.mean(data)) + """</td>
            </tr>
            <tr>
                <th>Mean Signal SD (unmasked)</th>
                <td>""" + str(np.std(data)) + """</td>
            </tr>
            <tr>
                <th>Mean Signal (masked)</th>
                <td>""" + str(np.mean(data[mask == 1])) + """</td>
            </tr>
            <tr>
                <th>Mean Signal SD (masked)</th>
                <td>""" + str(np.std(data[mask == 1])) + """</td>
            </tr>
        </table>"""
    html_output.write(mean_signal_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add Mask to html file
    html_output.write("<div id=\"Masks\"> <h1 class=\"sub-report-title\">Masks</h1>")
    masktext = """<p>Voxels included in the masks and used for the quality check are highlighted in blue,
                voxels used for SNR calculation are highlighted in green. 
                See <a href = "https://drmichaellindner.github.io/pyfMRIqc/#masks">here</a>
                for guidelines about how to use this plot.</p>"""
    html_output.write(masktext)
    html_output.write("""<img src="pyfMRIqc_MASK_""" + fname + """.png" alt="mask image" class="center">""")
    html_output.write("<p></p>")
    html_output.write("<h2> Mask summary</h2>")
    mean_table = """
        <table style="width:50%">
            <tr>
                <th>Total Number of Voxels</th>
                <td>""" + str(nrvoxel) + """</td>
            </tr>
            <tr>
                <th>Mask Threshold Value</th>
                <td>""" + str(mask_thresh) + """</td>
            </tr>
            <tr>
                <th>Number of Masked Voxels</th>
                <td>""" + str(np.sum(mask)) + """</td>
            </tr>
        </table>"""
    html_output.write(mean_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    html_output.write("""<div id="Variance"> <h1>Variance of voxel intensity</h1>""")
    html_output.write("<p>For better visualization the image is thresholded at max " + str(varthresh) +
                      " to minimize the scaling effect of large outliers.")
    vartext = """ See <a href = "https://drmichaellindner.github.io/pyfMRIqc/#variance-of-voxel-intensity">
                here</a> for guidelines about how to use this plot.</p>"""
    html_output.write(vartext)
    head, vfilename = os.path.split(varfilename)
    html_output.write("""<img src=""" + vfilename[:-4] + """_thr""" + str(varthresh) + vfilename[-4:] +
                      """ alt="Signal variance from functional image" class="center">""")

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add SNR data to html file
    html_output.write("""<div id="SNR"> <h1>Signal to noise ratio (SNR)</h1>""")
    snrtext = """<p>See <a href = "https://drmichaellindner.github.io/pyfMRIqc/#signal-to-noise-ratio">
                       here</a> for more info.</p><p></p>"""
    html_output.write(snrtext)

    allsclicesnrtext = ""
    for ii in range(len(snrvec)):
        s = str(snrvec[ii])
        s = s.replace('[', '')
        s = s.replace(']', '')
        s = s.replace(' ', '')
        allsclicesnrtext += s + "\n"

    snr_table = """
        <table style="width:70%">
            <tr>
                <th>Mask Threshold Value</th>
                <td>""" + str(mask_thresh) + """</td>
            </tr>
            <tr>
                <th>Number of Masked Voxels</th>
                <td>""" + str(np.sum(mask)) + """</td>
            </tr>
            <tr>
                <th>SNR Threshold</th>
                <td>""" + str(snr_voxel_percentage) + """</td>
            </tr>
            <tr>
                <th>Number of voxels below threshold for SNR</th>
                <td>""" + str(snrvoxel_n) + """</td>
            </tr>
            <tr>
                <th>Mean values of voxel for SNR</th>
                <td>""" + str(snrvoxel_mean) + """</td>
            </tr>
            <tr>
                <th>STD of voxel for SNR</th>
                <td>""" + str(snrvoxel_std) + """</td>
            </tr>
            <tr>
                <th>Value range of voxel for SNR</th>
                <td>""" + str(snrvoxel_min) + " - " + str(snrvoxel_max) + """</td>
            </tr>
                <tr>
                <th>Mean voxel SNR</th>
                <td>""" + str(np.nanmean(snrvec)) + """</td>
            </tr>
            <tr>
                <th>Min Slice SNR</th>
                <td>""" + str(np.nanmin(snrvec)) + """</td>
            </tr>
            <tr>
                <th>Max Slice SNR</th>
                <td>""" + str(np.nanmax(snrvec)) + """</td>
            </tr>
            <tr>
                <th>ALL Slice SNR</th>
                <td>""" + str(allsclicesnrtext) + """</td>
            </tr>
        </table>"""
    html_output.write(snr_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add Mask to html file
    html_output.write(
        "<div id=\"SUMDIFF\"> <h1 class=\"sub-report-title\">Sum of squared scaled difference over time</h1>")
    html_output.write(
        "<p>For better visualization the image is Mean squared scaled difference (SSD) is thresholded at max: " + str(mssdthresh) + ".")
    ssdtext = """This image shows the sum of the difference in voxel intensity between adjacent volumes for each 
                    voxel across all adjacent volumes acquired during scanning. See <a href = 
                    "https://drmichaellindner.github.io/pyfMRIqc/#sum-of-squared-scaled-differences-over-time">
                    here</a> for more info.</p>"""
    html_output.write(ssdtext)
    html_output.write("""<img src="pyfMRIqc_SUM_SQUARED_DIFF_SCALED_""" + fname +
                      """.png" alt="sum of squared scaled difference image" class="center">""")

    html_output.write("</div><br><hr><br>")  # horizontal line

    # Add motion parameter to html file
    html_output.write("""<div id="Motion"> <h1>Motion parameter summary</h1>""")

    if rp_fname is not None:
        motion_table = """
                <table style="width:70%">
                    <tr>
                        <th>Mean absolute Movement</th>
                        <td>""" + str(np.mean(rm)) + """</td>
                    </tr>
                    <tr>
                        <th>Max absolute Movement</th>
                        <td>""" + str(np.max(rm)) + """</td>
                    </tr>
                    <tr>
                        <th>Mean relative Movement</th>
                        <td>""" + str(np.mean(relrm)) + """</td>
                    </tr>
                    <tr>
                        <th>Max relative Movement</th>
                        <td>""" + str(np.max(relrm)) + """</td>
                    </tr>
                    <tr>
                        <th>Relative movements (>0.1mm)</th>
                        <td>""" + str(len(nrrm01)) + """</td>
                    </tr>
                    <tr>
                        <th><font color=#ffa500>Relative movements (>0.5mm)</font></th>
                        <td><font color=#ffa500>""" + str(len(nrrm05)) + """</font></td>
                    </tr>
                    <tr>
                        <th><font color="red">Relative movements (>voxelsize)</font></th>
                        <td><font color="red">""" + str(len(nrrmv)) + """</font></td>
                    </tr>
                </table>"""
    else:
        motion_table = """<p>No motion parameter file was provided as input.</p>"""

    html_output.write(motion_table)

    html_output.write("</div><br><hr><br>")  # horizontal line

    # get time string
    t = datetime.now()
    strg = t.strftime('%Y/%m/%d %H:%M:%S')

    html_output.write("""<div id="About"> <h1>About</h1>""")
    about_text = """
    <br>
    <ul>
		<li>Date of quality check: """ + str(strg) + """</li>
        <li>pyfMRIqc code: 
        <a href="https://github.com/DrMichaelLindner/pyfMRIqc">https://github.com/DrMichaelLindner/pyfMRIqc</a></li>
    </ul>
    <br>
    <p><font size="4"><b>Thank you for using pyfMRIqc.py!</b></font></p>
    <p><b>AUTHORS:</b><br>
        Michael Lindner and Brendan Williams<br>
        University of Reading, 2019<br>
        School of Psychology and Clinical Language Sciences<br>
        Center for Integrative Neuroscience and Neurodynamics
        </p>
    """
    html_output.write(about_text)
    html_output.write("</div><br><hr><br>")  # horizontal line

    # close html files
    html_output.write("</body></html>")
    html_output.close()

    print("DONE!")

def nii2image(img3d, cond, png_filename):

    # function from fmripyqc... 
    
    # nifti info
    slice_x, slice_y, slice_z = img3d.shape[0], img3d.shape[1], img3d.shape[2]
    matdim = np.ceil(np.sqrt(slice_z))
    img_x, img_y = int(matdim * slice_x), int(matdim * slice_y)

    # create image for plot
    image = np.zeros((img_y, img_x)) # flip x and y
    cx, cy = 0, 0
    for ii in reversed(range(slice_z)):
        xy = img3d[:, :, ii]
        image[cy:cy + slice_y, cx:cx + slice_x] = np.rot90(xy) 
        cx += slice_x
        if cx >= img_x:
            cx = 0
            cy += slice_y

    # specific details for diff. plots
    if cond == 'Variance':
        h = np.histogram(image, bins=1000)
        thr = h[1][max(np.argwhere(h[0] > 400))]
        vmax = np.round(thr[0], decimals=0)
        png_filename = f"{png_filename[:-4]}_thr{str(vmax)}{png_filename[-4:]}"
    elif cond == 'DIFF':
        hs, be = np.histogram(image, bins=50, density=True)
        vmax = np.round(be[5], decimals=0)
    else:
        vmax = img3d.max()

    sizes  = np.shape(image)
    height = float(sizes[0])
    width  = float(sizes[1])

    # plot
    fig = plt.figure()
    fig.set_size_inches(width / height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    ax.imshow(image, vmax=vmax, cmap='gray')
    if cond not in ['Mask', 'SNR']:
        plt.savefig(png_filename, dpi=height)
        plt.close()

    if cond == 'Variance': return vmax
    elif cond == 'DIFF':   return vmax
    else:                  return image

def get_roi_snr(snr_image, mask, **kwargs):

    # different options for mask
    if isinstance(mask, str):
        if mask in ['whole-brain', 'gm']:
            masker = NiftiMasker(mask_strategy=f'{mask}-template', **kwargs)
        else:
            masker = NiftiMasker(mask_img=mask, **kwargs)
    else: # should be a tuple of coordinates - doesnt return all voxels, just average
        masker = NiftiSpheresMasker([mask], **kwargs) 

    # mask the snr image 
    # maybe also return a nifti with just the snr? or plot?
    return masker.fit_transform(snr_image).T # voxels

def run_roi_snr(sub_id, img_dir, mask_dir):

    sub_images = [snr for snr in glob.glob(f'{img_dir}/*snr*.nii') if f'/{sub_id}_' in snr]
    out_dir = img_dir

    print(f'Running ROI SNR analysis for {sub_images}...')

    # check if summary already exists
    sub_fname = f'{out_dir}/{sub_id}_roi_snr_summary.xlsx'
    if os.path.exists(sub_fname):
        print(f'{sub_id}: already processed, skipping...')
    else: 

        sub_df = pd.DataFrame(columns=['sub_id', 'roi', 'snr_type', 'snr'])
        for mask in glob.glob(f'{mask_dir}/*.nii'): 
            
            # simplify mask name
            mask_nii  = mask.split('/')[-1]
            region    = ('_').join(mask_nii.split('_')[1:3])
            threshold = mask_nii.split('-thr')[-1].replace('-1mm.nii', '')
            mask_name = f'{region}_thr{threshold}'
            print(f'{sub_id}: getting SNR for {mask_name}', end='\r')

            # calculate snr for each image
            for snr_image in sub_images:
                print(f'{sub_id}: getting SNR for {mask_name} in {snr_image}', end='\r')
                snr_type = snr_image.split('/')[-1].split('_')[1]
                try: 
                    snr = np.mean(get_roi_snr(snr_image, mask))
                except:
                    snr = np.nan
                    print(f'Error: SNR for {sub_id} {mask_name} {snr_type}')
                sub_df.loc[len(sub_df), :] = [sub_id, mask_name, snr_type, snr]
        sub_df.to_excel(sub_fname, index=False)