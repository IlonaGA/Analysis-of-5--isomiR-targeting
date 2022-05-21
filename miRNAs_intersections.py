#%%
# This script was created to find intersections between given miRNAs

#%%
import pandas as pd 
import numpy as np
import os

# %%
miRNA_list = ['hsa-miR-93-5p|0'] #List example
dir = '/path/to/files/activity_analysis' #Path to results obtained by running 2_activity_analysis.py 

files = [file for file in os.listdir(dir) if file.split('_')[1] == 'tumor.tsv']

df = [pd.read_csv(dir + '/' + file, sep='\t', index_col=0)[['activity']]
            .rename(columns={'activity': file.split('_')[0]}) for file in files]
total_activity_df = df[0].loc[miRNA_list]

for i in range(1, len(df)):
    total_activity_df = total_activity_df.join(df[i], how='outer')
total_activity_df = total_activity_df.loc[miRNA_list]
total_activity_df.fillna(0, inplace=True)
total_activity_df = total_activity_df > 0 
total_activity_df = total_activity_df[np.sum(total_activity_df, axis=1) >= 2]
total_activity_df = total_activity_df.loc[:, np.sum(total_activity_df, axis=0) == total_activity_df.shape[0]]
projects = list(total_activity_df.columns)
print(projects)

# %%
target_dir = '/path/to/files/corr_analysis' #Path to results obtained by running 1_create_raw_and_summary.py 
target_files = [target_file for target_file in os.listdir(target_dir) 
                if target_file.split('_')[2] == 'raw.tsv' and 
                   target_file.split('_')[1] == 'tumor' and
                   target_file.split('_')[0] in projects]

for target_file in target_files:
    miRNA1_set = pd.read_csv(target_dir + '/' + target_file, sep='\t', index_col=0)[['gene', 'corr']].loc[miRNA_list[0]]
    miRNA1_set = miRNA1_set[miRNA1_set['corr'] <= -0.3]
    miRNA1_set = set(list(miRNA1_set['gene']))
    for i in range (1, len(miRNA_list)):
        target_df = pd.read_csv(target_dir + '/' + target_file, sep='\t', index_col=0)[['gene', 'corr']].loc[miRNA_list[i]]
        target_df = target_df[target_df['corr'] <= -0.3]
        target_set = set(list(target_df['gene']))
        intersection = miRNA1_set.intersection(target_set)
        print(target_file, intersection)
