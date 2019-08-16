# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/07/19
content:    Preprocess brain data from Darmanis et al. 2015/2017.
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    fn_brain_atlas_meta = '../data/brain_atlas/GSE67835_family.soft'
    fns_brain_atlas_counts = glob.glob('../data/brain_atlas/GSM*.csv')

    print('Parse brain atlas metadata')
    cells = {'name': [], 'cellType': []}
    with open(fn_brain_atlas_meta, 'rt') as f:
        for line in f:
            if line.startswith('^SAMPLE ='):
                name = line.rstrip('\n').split()[-1]
                cells['name'].append(name)
            elif line.startswith('!Sample_characteristics_ch1 = cell type:'):
                ctype = line.rstrip('\n').split()[-1]
                cells['cellType'].append(ctype)
    cells = pd.DataFrame(cells).set_index('name')

    print('A bunch of cell types are misannotated')
    cell_types = cells['cellType'].values
    cell_types[cell_types == 'astrocytes'] = 'Astrocyte'
    cell_types[cell_types == 'neurons'] = 'Neuron'
    cell_types[cell_types == 'endothelial'] = 'Vascular'
    cell_types[cell_types == 'oligodendrocytes'] = 'Oligodendrocyte'
    cells['cellType'] = cell_types

    print('Exclude hybrids, whatever they are')
    ind = cell_types != 'hybrid'
    cells = cells.loc[ind]

    print('Parse brain atlas counts')
    ncells = len(cells)
    fns_dict = {x.split('/')[-1].split('_')[0]: x for x in fns_brain_atlas_counts}
    counts = None
    for ic, cn in enumerate(cells.index):
        count = pd.read_csv(fns_dict[cn], sep='\t', index_col=0, squeeze=True)
        if counts is None:
            ngenes = count.shape[0]
            genes = count.index.str.rstrip(' ')
            counts = np.zeros((ngenes, ncells), np.float32)
        counts[:, ic] = count.values
    counts = pd.DataFrame(counts, index=genes, columns=cells.index)

    print('Store brain atlas to file')
    loompy.create(
        '../data/brain_atlas/dataset.loom',
        layers={'': counts.values},
        row_attrs={'GeneName': counts.index.values},
        col_attrs={'CellID': counts.columns.values, 'cellType': cells['cellType'].values},
        )

    print('Parse GBM metadata')
    fn_gbm_meta = '../data/glioblastoma/GSE84465_family.soft'
    fn_gbm_counts = '../data/glioblastoma/GSE84465_GBM_All_data.csv'

    print('Parse GBM metadata')
    cells_gbm = {'name': [], 'cellType': []}
    with open(fn_gbm_meta, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_description = 1001'):
                name = line.rstrip('\n').split()[-1]
                cells_gbm['name'].append(name)
            elif line.startswith('!Sample_characteristics_ch1 = cell type:'):
                ctype = line.rstrip('\n').split()[-1]
                cells_gbm['cellType'].append(ctype)
    cells_gbm = pd.DataFrame.from_dict(cells_gbm).set_index('name')
    ncells_gbm = len(cells_gbm)

    print('A bunch of cell types are misannotated')
    cell_types = cells_gbm['cellType'].values
    cell_types[cell_types == 'Astocyte'] = 'Astrocyte'
    cell_types[cell_types == 'Neoplastic'] = 'Neoplastic1'
    cell_types[cell_types == 'cell'] = 'Neoplastic2'
    cells_gbm['cellType'] = cell_types

    print('Parse GBM counts')
    counts_gbm = pd.read_csv(fn_gbm_counts, sep=' ', index_col=0).astype(np.float32)

    print('Store GBM to file')
    loompy.create(
        '../data/glioblastoma/dataset.loom',
        layers={'': counts_gbm.values},
        row_attrs={'GeneName': counts_gbm.index.values},
        col_attrs={'CellID': counts_gbm.columns.values, 'cellType': cells_gbm['cellType'].values},
        )

    print('Find intersection of genes')
    genes_int = np.intersect1d(genes, counts_gbm.index.values)

    print('Merge datasets without averaging')
    cells_merge = pd.Index(cells.index.tolist() + cells_gbm.index.tolist(), name='CellID')
    cells_dataset = ['atlas'] * ncells + ['new'] * ncells_gbm
    cells_dataset = np.array(cells_dataset)
    cell_types = cells['cellType'].values.tolist() + cells_gbm['cellType'].values.tolist()
    cell_types = np.array(cell_types)
    counts_merge = np.empty((len(genes_int), len(cells_merge)), dtype=np.float32)
    counts_merge[:, :len(cells)] = counts.loc[genes_int].values
    counts_merge[:, len(cells):] = counts_gbm.loc[genes_int].values
    counts_merge = pd.DataFrame(counts_merge, index=genes_int, columns=cells_merge)

    print('Filter and normalize joint dataset')
    # Last 4 features are trash
    counts_merge = counts_merge.iloc[:-4]
    cov = counts_merge.sum(axis=0)
    counts_merge *= 1e6 / cov

    print('Store both brains to file, no averaging')
    loompy.create(
        '../data/both_brain/dataset_full.loom',
        layers={'': counts_merge.values},
        row_attrs={'GeneName': counts_merge.index.values},
        col_attrs={
            'CellID': counts_merge.columns.values,
            'coverage': cov.values,
            'cellType': cell_types,
            'dataset': cells_dataset,
            },
        )

    print('Average atlas')
    counts = counts.T
    counts['cellType'] = cells['cellType']
    avgs = counts.groupby('cellType').mean().T

    print('Merge datasets with averages')
    newindex = avgs.columns.tolist() + cells_gbm.index.tolist()
    counts_merge_avg = counts_gbm.loc[genes_int].copy()
    cells_merge_avg = cells_gbm.copy()
    for col in avgs:
        counts_merge_avg[col] = avgs.loc[genes_int, col]
        cells_merge_avg.loc[col] = {'cellType': col}
    counts_merge_avg = counts_merge_avg[newindex]
    cells_merge_avg = cells_merge_avg.loc[newindex]

    print('Filter and normalize joint dataset')
    # Last 4 features are trash
    counts_merge_avg = counts_merge_avg.iloc[:-4]
    cov = counts_merge_avg.sum(axis=0)
    counts_merge_avg *= 1e6 / cov

    print('Store both brains to file, atlas averaged')
    loompy.create(
        '../data/both_brain/dataset.loom',
        layers={'': counts_merge_avg.values},
        row_attrs={'GeneName': counts_merge_avg.index.values},
        col_attrs={
            'CellID': counts_merge_avg.columns.values,
            'coverage': cov.values,
            'cellType': cells_merge_avg['cellType'].values,
            },
        )
