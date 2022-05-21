import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy import stats
import sys
import os

sys.path.insert(1, '/path/to/directory/') #Path to working directory

from tqdm import tqdm
import core
from core.correlation_utils import spearmanr_thr_proba


def molecule_activity(name, interaction, thr=-0.3):
    indexes = interaction['isomiR'] == name
    corr = np.sum(interaction['corr'][indexes] <= thr)
    return corr 


def count_activity(interaction, thr=-0.3):
    activity_table = []

    for name in pd.unique(interaction['isomiR']):
        corr = molecule_activity(name, interaction, thr)
        mirna_median = np.median(data_df.loc[name])
        number_of_targets = interactions_dict[name]
        activity_table.append([name, mirna_median, number_of_targets, corr])

    activity_table = pd.DataFrame(activity_table, columns=['isomiR', 'expression_median', 'predicted_targets', 'sign_corr'])
    return activity_table


TCGA_project = sys.argv[1]
interaction_path = sys.argv[2] #path should not contain '/' in the end

threshold = -0.3

data_df = pd.read_csv(f'/path/to/files/{TCGA_project}_RPM.tsv', sep='\t', index_col=0) #Path to RPM files
interaction_df = pd.read_csv(f'{interaction_path}/path/to/files/{TCGA_project}.tsv', sep='\t') #Path to files with tagrets
interaction_df.columns = ['isomiR', 'gene']

interactions_dict = {}
for miRNA in interaction_df['isomiR'].unique():
    interactions_dict[miRNA] = np.sum(interaction_df['isomiR'] == miRNA)

interaction_df['corr'] = core.spearmanr(data_df, interaction_df['isomiR'], interaction_df['gene'])
interaction_df = interaction_df.sort_values(by='corr')

interaction_df['proba'] = spearmanr_thr_proba(threshold, interaction_df['corr'], data_df.shape[1])

activity_table = count_activity(interaction_df)

interaction_df['sample size'] = data_df.shape[1]

expression_medians = pd.DataFrame(data_df.median(axis=1))
expression_medians.columns = ["isomiR_expression_median"]
interaction_df = interaction_df.set_index("isomiR").join(expression_medians).reset_index()
interaction_df = interaction_df.rename(columns={"index": "isomiR"})
expression_medians.columns = ["gene_expression_median"]
interaction_df = interaction_df.set_index("gene").join(expression_medians).reset_index()
interaction_df = interaction_df.rename(columns={"index": "gene"})
interaction_df = interaction_df[['isomiR', 'gene', 'isomiR_expression_median', 'gene_expression_median', 'corr', 'sample size', 'proba']]

print('END ' + TCGA_project)
interaction_df.to_csv(f'{interaction_path}/corr_analysis/{TCGA_project}_raw.tsv', index=None, sep='\t')
activity_table.to_csv(f'{interaction_path}/corr_analysis/{TCGA_project}_summary.tsv', index=None, sep='\t')
