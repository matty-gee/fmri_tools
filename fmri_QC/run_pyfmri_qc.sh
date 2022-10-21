#!/bin/bash

# https://openresearchsoftware.metajnl.com/articles/10.5334/jors.280/
# https://drmichaellindner.github.io/pyfMRIqc/
# python pyfMRIqc.py \
#           -n <your_raw_functional_file.nii> \
#           -s % of voxels outside of mask w/ lowest mean time course values used to compute SNR \
#           -k <your_brain_mask_file.nii> \
#           -o <your_output_path> \
#           -m <your_motion_file> \
#           -x

# mask img
mask_img=./masks/EPI_resampled_thr20.nii
s=10
t=1000
for SUBDIR in subs/*; 
do
    if [ $SUBDIR != 'subs/@eaDir' ]; then
        # func img
        func_img=$(find $SUBDIR/func -type f -iname "2*.nii") # wants raw image
    
        # motion parameters
        rp_txt=$(find $SUBDIR/func -type f -iname "rp*.txt")

        # output dir 
        mkdir -p $SUBDIR/pyfmri_qc # make only if doesnt exist
        outdir=$SUBDIR/pyfmri_qc

        # # run the command 
        if [ -f $outdir/*.html ]; then
            echo "subject has been run"
        else
            echo "running for" $func_img
            python ../../../../Code/fMRI/qc/pyfMRIqc-master/pyfMRIqc.py \
                -n $func_img \
                -s $s \
                -t $t \ # I need a mask image w/ same shape... for now using this intensity threshold
                -o $outdir \
                -m $rp_txt 
        fi
    fi
done