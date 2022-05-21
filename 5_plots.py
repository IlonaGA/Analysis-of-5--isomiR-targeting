#%%
import sys
import os
import pickle

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

from scipy.stats import *

import matplotlib.colors as mcol
import matplotlib as mpl

import matplotlib.gridspec as gridspec
import SeabornFig2Grid as sfg


#%%
def rc(path, i=0, h="infer"):
    return pd.read_csv(path, index_col=i, header=h)


def rt(path, i=0, h="infer"):
    return pd.read_csv(path, sep="\t", index_col=i, header=h)



#%%
#creating DF
databases = ['miRDB', 'TargetScan', 'miRDB_and_TargetScan']

for database in databases:
    files = []
    dotplot_df = []

    for file in os.listdir(f'/huge/steve/isomiR_research/{database}/corr_analysis/TMM-RPM/'):
        files.append(file.split('_')[0])
    files = set(files)

    for project in files:
        data = pd.read_csv(f'/huge/steve/isomiR_research/{database}/corr_analysis/TMM-RPM/{project}_tumor_summary.tsv', sep='\t', index_col=0)
        highly_expressed = pd.read_csv(f'/huge/steve/isomiR_research/expression_data/highly_expressed_isomiRs/{project}_tumor.tsv', sep='\t', index_col=0)

        data = data.loc[set(data.index) & set(highly_expressed.index)]

        spearman = spearmanr(data['sign_corr'], data['predicted_targets'])
        dotplot_df.append([project, spearman[0], len(highly_expressed)])
    
    dotplot_df = pd.DataFrame(dotplot_df, columns=['project', 'corr', 'highly_expressed'])
    dotplot_df.to_csv(f'{database}_df.tsv', sep='\t', index=None)

#%%
mpl.rcParams["font.sans-serif"] = "Arial"

def dot_plot(in_fname, out_fname):
    df = rt(in_fname, i=None)
    
    df = df.sort_values('corr')
    df = df.reset_index().reset_index()

    plt.figure(figsize=(6.5, 10))

    cmap = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","b"])

    ax = sns.scatterplot(
        x='corr', y='level_0', hue='highly_expressed', data=df,
        palette=cmap,
        #hue_norm=(df['highly_expressed'].min(), 0.05),
        legend=False
    )

    norm = plt.Normalize(df['highly_expressed'].min(), df['highly_expressed'].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    ax.figure.colorbar(sm, label='Number of highly expressed miRNAs')

    plt.xlabel('Spearman correlation')
    plt.ylabel('')
    
    plt.xlim([df['corr'].min()*0.9, df['corr'].max()*1.1])
    plt.yticks(ticks=df['level_0'], labels=df['project'], size=6)

    plt.savefig(out_fname)


#%%
#SCATTERPLOT
databases = ['miRDB_and_TargetScan']

for database in databases:
    files = ['TCGA-CHOL', 'TCGA-PCPG']
    
    #project 1
    project = files[0]
    data = pd.read_csv(f'/huge/steve/isomiR_research/{database}/corr_analysis/TMM-RPM/{project}_tumor_summary.tsv', sep='\t', index_col=0)
    highly_expressed = pd.read_csv(f'/huge/steve/isomiR_research/expression_data/highly_expressed_isomiRs/{project}_tumor.tsv', sep='\t', index_col=0)

    data = data.loc[set(data.index) & set(highly_expressed.index)]

    spearman = spearmanr(data['sign_corr'], data['predicted_targets'])
    spearman_str = 'Spearman correlation = {:.3f}'.format(spearman[0])
    data.rename(columns={'predicted_targets':'Number of predicted targets', 'sign_corr': f'Number of anti-correlated targets, {project} \n{spearman_str}'}, inplace=True)
    p0 = sns.jointplot(data=data, x=f'Number of anti-correlated targets, {project} \n{spearman_str}', y='Number of predicted targets', marginal_kws={'bins':10})
    p0.fig.suptitle('Spearman = {:2f} p_value = {:2f}, {}'.format(spearman[0], spearman[1], project))
    p0.fig.subplots_adjust(top=0.95)

    #project 2
    project = files[1]
    data = pd.read_csv(f'/huge/steve/isomiR_research/{database}/corr_analysis/TMM-RPM/{project}_tumor_summary.tsv', sep='\t', index_col=0)
    highly_expressed = pd.read_csv(f'/huge/steve/isomiR_research/expression_data/highly_expressed_isomiRs/{project}_tumor.tsv', sep='\t', index_col=0)

    data = data.loc[set(data.index) & set(highly_expressed.index)]

    spearman = spearmanr(data['sign_corr'], data['predicted_targets'])
    spearman_str = 'Spearman correlation = {:.3f}'.format(spearman[0])
    data.rename(columns={'predicted_targets':'Number of predicted targets', 'sign_corr': f'Number of anti-correlated targets, {project} \n{spearman_str}'}, inplace=True)
    p1 = sns.jointplot(data=data, x=f'Number of anti-correlated targets, {project} \n{spearman_str}', y='Number of predicted targets', marginal_kws={'bins':10})
    p1.fig.suptitle('Spearman = {:2f} p_value = {:2f}, {}'.format(spearman[0], spearman[1], project))
    p1.fig.subplots_adjust(top=0.95)

    fig = plt.figure(figsize=(6.5, 10))
    gs = gridspec.GridSpec(2, 1)
    
    img0 = sfg.SeabornFig2Grid(p0, fig, gs[0]) 
    img1 = sfg.SeabornFig2Grid(p1, fig, gs[1]) 

    gs.tight_layout(fig)

    #dot_plot('miRDB_and_TargetScan_df.tsv', 'miRDB_and_TargetScan_dotplot.pdf', 0.5)
    #plt.show()
    plt.savefig('jointplots.png')
    dot_plot('miRDB_and_TargetScan_df.tsv', 'miRDB_and_TargetScan_dotplot.png')
    ! convert jointplots.png miRDB_and_TargetScan_dotplot.png +append res.png

# %%
data_a = pd.DataFrame(np.random.rand(200).reshape(100, 2), columns=['x', 'y'])
data_a['source'] = 'a'
data_b = pd.DataFrame(np.random.rand(200).reshape(100, 2), columns=['x', 'y'])
data_b['source'] = 'b'

data = pd.concat([data_a, data_b])
data.reset_index(inplace=True)

g0 = sns.jointplot(data=data[data['source'] == 'a'], x='x', y='y')
g1 = sns.jointplot(data=data[data['source'] == 'b'], x='x', y='y')


fig = plt.figure(figsize=(5,13))
gs = gridspec.GridSpec(2, 1)

mg0 = sfg.SeabornFig2Grid(g0, fig, gs[0])
mg1 = sfg.SeabornFig2Grid(g1, fig, gs[1])

gs.tight_layout(fig)

plt.show()


