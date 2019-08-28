# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/08/19
content:    Preprocess human liver by Aizarani et al. 2019
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    fn_atlas_meta = '../data_raw/Aizaran_2019/GSE124395_clusterpartition.tsv.gz'
    fn_atlas_counts = '../data_raw/Aizaran_2019/GSE124395_Normalhumanlivercellatlasdata.tsv.gz'

    print('Parse atlas metadata')
    meta = pd.read_csv(fn_atlas_meta, sep=' ', compression='gzip', index_col=0)

    print('Parse atlas counts')
    counts = pd.read_csv(fn_atlas_counts, sep='\t', compression='gzip', index_col=0)

    # Some cells have no metadata
    counts = counts.loc[:, meta.index].astype(np.float32)

    # Reconstruction from the horrible t-SNE plus marker genes
    cell_typed = {
        'Hepatocyte': [11, 14, 17, 30],
        'Liver sinusoidal endothelial cell': [9, 13, 20],
        'Macrovascular endothelial cells': [10, 29],
        'Epithelial cell': [4, 7, 24, 39],
        'B cell': [8, 22, 38],
        'Kupffer cell': [2, 6, 23, 25, 31],
        'T/NK cell': [1, 3, 5, 12, 18, 28],
        'Stellate cell/Myofibroblast': [33],
        }
    meta['CellType'] = 'Unknown'
    for ct, clns in cell_typed.items():
        meta.loc[meta.iloc[:, 0].isin(clns), 'CellType'] = ct

    print('Exclude unknown')
    ind = meta['CellType'] != 'Unknown'
    counts = counts.loc[:, ind]
    meta = meta.loc[ind]

    print('Store atlas to file')
    os.makedirs('../data_full/Aizaran_2019')
    loompy.create(
        '../data_full/Aizaran_2019/dataset.loom',
        layers={'': counts.values},
        row_attrs={'GeneName': counts.index.values},
        col_attrs={
            'CellID': counts.columns.values,
            'cellType': meta['CellType'].values},
        )
