'''
'''

import sys
import pandas as pd
import numpy as np
from generic import read_excel

sys.path.insert(0, f'{code_dir}/fMRI/mni_to_region')
from find_mni import *

shen_df = read_excel('/Users/matty_gee/Desktop/Shen_parcellations/Shen_368/Shen_368_labels.xlsx')
mni_coords = shen_df[['MNI_X','MNI_Y','MNI_Z']]

# one approach...
[one_line, table] = find_structure(mni_coords[1:])
region_names = [t[0].replace(' ', '_') for t in table]
region_names = pd.DataFrame(['Background'] + region_names, columns=['Region_name']).astype(str)
shen_df2 = pd.concat([shen_df, region_names], 1)
for reg in np.unique(shen_df2['Region_name']): 
    mask = shen_df2['Region_name'] == reg
    shen_df2.loc[mask,'Region_name'] = [reg + '_' + str(i) for i in np.arange(1, np.sum(mask)+1)]
shen_df2.to_excel('/Users/matty_gee/Desktop/Shen_parcellations/Shen_368/Shen_368_labels2.xlsx', index=False)