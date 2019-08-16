# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/19
content:    Preprocess B cells from Croote et al. 2019, only half the
            individuals to reduce batch effects
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    print('Parse B cell metadata')
    fn_pbmc_meta = '../data/Bcells_croote/croote_bcells_naivemem_PB_meta.csv'
    cells = pd.read_csv(fn_pbmc_meta, sep=',', index_col=0)
    cells['cellType'] = ''
    cells.loc[cells['cluster'] == 0, 'cellType'] = 'B cell, naive_memory'
    cells.loc[cells['cluster'] == 1, 'cellType'] = 'B cell, plasmablast'

    print('Parse count data')
    fn_pbmc_counts = '../data/Bcells_croote/croote_bcells_naivemem_PB_raw_cnts.csv.gz'
    counts = pd.read_csv(
        fn_pbmc_counts,
        sep=',',
        compression='gzip',
        index_col=0,
        )
    counts = counts.astype(np.float32)

    print('Store B cell data to file')
    loompy.create(
        '../data/Bcells_croote/dataset.loom',
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
            'Number of cells': counts.shape[1],
            }
        )
