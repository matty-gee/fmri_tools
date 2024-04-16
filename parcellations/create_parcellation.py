import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import plotting
from six.moves import cPickle as pickle 
def pickle_file(file_, filename_):
    with open(filename_, 'wb') as f:
        pickle.dump(file_, f)
    f.close()
def load_pickle(filename_):
    with open(filename_, 'rb') as f:
        ret_file = pickle.load(f)
    return ret_file

#####################################################################
## make parcellation example...
#####################################################################

parcellation_dir = f'{user}/Desktop/Shen_parcellations/Shen_368'
labels_df = pd.read_excel(f'{parcellation_dir}/Shen_368_labels.xlsx')
labels_df = labels_df.sort_values(by='node_index')
display(labels_df.head(5))

# show whole parcellation
fname = f'{parcellation_dir}/Shen_1mm_368_parcellation.nii.gz'
parcellation_nii = nib.load(fname)
plotting.plot_roi(parcellation_nii)
plt.show()

# check specific BAs to check that labels match
parcellation_data = parcellation_nii.get_fdata()
node_nii = nib.Nifti1Image((parcellation_data == 0) * 1, parcellation_nii.affine)
plotting.plot_roi(node_nii, cut_coords=(-30,0,-36))
plt.show()

# pickle the parcellation
from sklearn.utils import Bunch
parcellation = Bunch(maps=parcellation_nii, labels=list(labels_df['BA']))
out_fname = f'{parcellation_dir}/Shen_1mm_368_parcellation.pkl'
pickle_file(parcellation, out_fname)

# reload pickle and plot again
parcellation_pkl = load_pickle(out_fname)
plotting.plot_roi(parcellation_pkl.maps)

#####################################################################