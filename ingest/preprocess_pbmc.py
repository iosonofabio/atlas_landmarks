# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/19
content:    Preprocess PBMCs from Zanini et al. 2018, only healthy individuals
'''
import os
import sys
import glob
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import loompy

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, CountsTable, FeatureSheet, SampleSheet


if __name__ == '__main__':

    fn_pbmc_meta = '../data/pbmc_zanini/samplesheet_10_10_unique_L1.tsv'

    print('Parse PBMC metadata')
    cells = pd.read_csv(fn_pbmc_meta, sep='\t', index_col='name')
    # Limit to uninfected individuals
    cells = cells.loc[cells['patient'].str.startswith('3-')]

    print('Parse count data')
    fn_pbmc_counts = '../data/pbmc_zanini/counts_10_10_unique_L1.tsv.gz'
    counts = CountsTable(pd.read_csv(
        fn_pbmc_counts,
        sep='\t',
        compression='gzip',
        index_col=0),
        )
    counts._normalized = 'counts_per_million'
    counts = counts.loc[:, cells.index].astype(np.float32)

    print('Parse gene metadata')
    fn_pbmc_meta_genes = '../data/pbmc_zanini/featuresheet_10_10_unique_L1.tsv'
    featuresheet = pd.read_csv(fn_pbmc_meta_genes, sep='\t', index_col=0)
    c = Counter(featuresheet['GeneName'])
    ind = []
    for fea in featuresheet.index:
        if c[featuresheet.at[fea, 'GeneName']] == 1:
            ind.append(fea)
    featuresheet = featuresheet.loc[ind]
    counts = counts.loc[ind]

    print('Reannotate cell types')
    ds = Dataset(
        samplesheet=cells,
        counts_table=counts,
        featuresheet=featuresheet,
        )
    ds.reindex(axis='features', column='GeneName', inplace=True)

    print('Filter low-quality cells')
    ds.samplesheet['coverage'] = ds.samplesheet['coverage'].astype(int)
    ds.query_samples_by_metadata('coverage >= 50000', inplace=True)

    print('Ignore HLA types and TCR and BCR variable regions')
    def fun(fn):
        if fn in ('HLA-A', 'HLA-B', 'HLA-C'):
            return False
        if fn.startswith('TRBV'):
            return False
        if fn.startswith('IGHV'):
            return False
        if fn.startswith('IGLV'):
            return False
        if fn.startswith('IGKV'):
            return False
        return True
    ds.counts = ds.counts.loc[list(filter(fun, ds.featurenames))]

    print('Cluster and embed')
    features = ds.feature_selection.overdispersed_within_groups(
        groupby='patient',
        n_features=500,
        inplace=False,
        )
    dsf = ds.query_features_by_name(features)
    dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')
    vs = dsc.dimensionality.tsne(perplexity=30)
    neighbors = dsc.graph.knn()
    edges = [tuple(y) for y in set(frozenset(x) for x in zip(neighbors.row, neighbors.col))]
    ds.samplesheet['cluster'] = [str(x) for x in dsc.cluster.leiden('samples', edges, resolution_parameter=0.005)]

    print('Calculate averages of genes by cluster')
    marker_genes = [
        'PTPRC', 'CD2', 'CD3E', 'GNLY', 'BATF3',
        'MS4A1', 'JCHAIN', 'IGHG1', 'TCL1A', 'PTPRS',
        'CD14', 'FCGR3A', 'CD1C', 'PPBP']
    clu = np.sort(ds.samplesheet['cluster'].unique())
    avg = pd.DataFrame(
        np.zeros((len(marker_genes), len(clu)), dtype=np.float32),
        index=marker_genes,
        columns=clu,
        )
    dsms = ds.query_features_by_name(marker_genes).split('cluster')
    for ui in clu:
        avg[ui] = dsms[ui].counts.mean(axis=1)
    avg1 = avg.copy()
    cell_types = defaultdict(list)
    cell_types['B cell'].append(avg1.loc['MS4A1'].idxmax())
    del avg1[cell_types['B cell'][-1]]
    cell_types['NK cell'].append(avg1.loc['GNLY'].idxmax())
    del avg1[cell_types['NK cell'][-1]]
    cell_types['pDC'].append(avg1.loc['PTPRS'].idxmax())
    del avg1[cell_types['pDC'][-1]]
    cell_types['nonclassical monocyte'].append(avg1.loc['FCGR3A'].idxmax())
    del avg1[cell_types['nonclassical monocyte'][-1]]
    # This one often makes two subclusters
    for i in range(2):
        cell_types['classical monocyte'].append(avg1.loc['CD14'].idxmax())
        del avg1[cell_types['classical monocyte'][-1]]
    cell_types['cDC'].append(avg1.loc['CD1C'].idxmax())
    del avg1[cell_types['cDC'][-1]]
    cell_types['T cell'].append(avg1.loc['CD3E'].idxmax())
    del avg1[cell_types['T cell'][-1]]
    cell_types['platelet'].append(avg1.loc['PPBP'].idxmax())
    del avg1[cell_types['platelet'][-1]]

    ds.samplesheet['cellType'] = ds.samplesheet['cluster']
    for ct, uis in cell_types.items():
        for ui in uis:
            ds.samplesheet.loc[ds.samplesheet['cluster'] == ui, 'cellType'] = ct

    print('Plot')
    marker_genes = [
        'PTPRC', 'CD2', 'CD3E', 'GNLY', 'BATF3',
        'MS4A1', 'JCHAIN', 'IGHG1', 'TCL1A', 'PTPRS',
        'CD14', 'FCGR3A', 'CD1C', 'PPBP']
    markers = marker_genes + ['cov+1', 'patient', 'cluster', 'cellType']
    ds.samplesheet['cov+1'] = ds.samplesheet['coverage'] + 1
    fig, axs = plt.subplots(3, 6, figsize=(13.5, 8), sharex=True, sharey=True)
    axs = axs.ravel()
    for gene, ax in zip(markers, axs):
        ax.set_title(gene)
        if gene in ('patient', 'cluster', 'cellType'):
            u = ds.samplesheet[gene].unique()
            lu = len(u)
            cmap = dict(zip(u, sns.color_palette('husl', n_colors=lu)))
            color_log = False
        else:
            u = None
            lu = None
            cmap = 'viridis'
            color_log = True
        ds.plot.scatter_reduced(
            vs,
            ax=ax,
            alpha=0.1,
            cmap=cmap,
            color_log=color_log,
            color_by=gene,
            s=12,
            )
        if gene in ('cluster', 'cellType'):
            for ui in u:
                vsc = vs.loc[ds.samplesheet[gene] == ui]
                xc, yc = vsc.values.mean(axis=0)
                ax.scatter([xc], [yc], s=50, facecolor='none', edgecolor=cmap[ui], marker='^', lw=4)
                ax.text(xc, yc, str(ui), fontsize=10, ha='center', va='bottom', clip_on=True)

    fig.tight_layout()
    plt.ion()
    plt.show()

    cells = ds.samplesheet[['cellType', 'patient']]
    counts = ds.counts

    print('Store PBMC data to file')
    loompy.create(
        '../data/pbmc_zanini/dataset.loom',
        layers={'': counts.values},
        row_attrs={'GeneName': counts.index.values},
        col_attrs={
            'CellID': counts.columns.values,
            'cellType': cells['cellType'].values,
            'individual': cells['patient'].values,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Species': 'Homo Sapiens',
            'Tissue': 'blood',
            'Age': 'adult',
            'Number of cells': ds.n_samples,
            }
        )
