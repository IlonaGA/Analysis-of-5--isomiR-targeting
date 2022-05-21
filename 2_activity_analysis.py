#%%
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
from scipy.stats import rankdata
from scipy.stats import hypergeom

# %%
TCGA_project = sys.argv[1]
folder = sys.argv[2]
miRDB = pd.read_csv(f'/path/to/files/{folder}/corr_analysis/{TCGA_project}.tsv', sep='\t') #Path to the directory with results obtained by running 1_create_raw_and_summary.py
neg_control = pd.read_csv(f'/path/to/files/negative_controls/corr_analysis/{TCGA_project}.tsv', sep='\t') #Path to the directory with results obtained by running 1_create_raw_and_summary.py
neg_control.head()

# %%
#filtering highly expressed isomiR
TCGA_project_for_filtration = TCGA_project.split('_')[0]
highly_expressed_isomiRs = pd.read_csv('../supp_data/highly_expressed_isomiRs.tsv', sep='\t') 
miRDB = miRDB.set_index('isomiR')
neg_control = neg_control.set_index('isomiR')
highly_expressed_isomiRs = highly_expressed_isomiRs[highly_expressed_isomiRs['project'] == TCGA_project_for_filtration]

#%%
miRDB_index = set(miRDB.index) & set(highly_expressed_isomiRs['isomiR'])
miRDB = miRDB.loc[miRDB_index]

neg_control_index = set(neg_control.index) & set(highly_expressed_isomiRs['isomiR'])
neg_control = neg_control.loc[neg_control_index]

miRDB = miRDB.reset_index()
neg_control = neg_control.reset_index()
neg_control.head()
#filtration ended

# %%
assert set(miRDB['isomiR']).issubset(set(neg_control['isomiR']))
neg_control.set_index('isomiR', inplace=True)
miRDB.set_index('isomiR', inplace=True)

#%%
activity_array_sign_corr = []

for isomir in miRDB.index:
    M = neg_control.loc[isomir, 'predicted_targets']
    n = neg_control.loc[isomir, 'sign_corr']
    N = miRDB.loc[isomir, 'predicted_targets']
    k = miRDB.loc[isomir, 'sign_corr']

    quantile = hypergeom(M, n, N).ppf(0.95)
    p_value = hypergeom(M, n, N).sf(k - 1)

    activity_array_sign_corr.append([isomir, 
                                     miRDB.loc[isomir, 'expression_median'],
                                     max(0, miRDB.loc[isomir, 'sign_corr'] - quantile),
                                     p_value])

TCGA_project_corr = pd.DataFrame(activity_array_sign_corr, columns=['isomiR', 'expression_median', 'activity', 'p_value'])
TCGA_project_corr['FDR'] = np.minimum(TCGA_project_corr['p_value'] * len(TCGA_project_corr) / rankdata(TCGA_project_corr['p_value']), 1)
TCGA_project_corr = TCGA_project_corr.sort_values('FDR')
TCGA_project_corr.to_csv(f'/path/to/files/{folder}/activity_analysis/{TCGA_project[:8]}_sign_corr.tsv', sep='\t', index=None)

