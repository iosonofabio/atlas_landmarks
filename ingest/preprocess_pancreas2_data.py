# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/07/19
content:    Preprocess pancreas data from Stitzel et al. 2016.
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    fn_pancreas_atlas_meta = '../data/pancreas_atlas2/GSE86469_family.soft'

    print('Parse brain atlas metadata (Stitzel et al. 2016)')
    cells = {'name': [], 'cellType': [], 'disease': []}
    with open(fn_pancreas_atlas_meta, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_title ='):
                name = line.rstrip('\n').split()[-1]
                cells['name'].append(name)
            elif line.startswith('!Sample_characteristics_ch1 = cell type:'):
                ctype = line.rstrip('\n').split(':')[-1][1:]
                cells['cellType'].append(ctype)
            elif line.startswith('!Sample_characteristics_ch1 = disease: '):
                dise = line.rstrip('\n').split(':')[-1][1:]
                cells['disease'].append(dise)
    cells = pd.DataFrame(cells).set_index('name')

    print('Exclude unsure cell types and diabetes patients')
    ind = cells['cellType'] != 'None/Other'
    ind &= cells['disease'] == 'Non-Diabetic'
    cells = cells.loc[ind]

    print('Download and parse pancreas atlas counts')
    ncells = len(cells)
    counts = pd.read_csv(
        '../data/pancreas_atlas2/GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz',
        sep=',',
        compression='gzip',
        index_col=0)
    counts = counts.loc[:, cells.index]

    print('Translate gene names')
    from collections import Counter
    tra = pd.read_csv(
        '../data/pancreas_atlas2/mart_export.tsv',
        sep='\t',
        index_col=0,
        squeeze=True,
        )
    c = Counter(tra.values)
    ind = [x for x in tra.index if c[tra.at[x]] == 1]
    tra = tra.loc[ind]
    ind = counts.index.isin(tra.index)
    counts = counts.loc[ind]
    idx = tra.loc[counts.index]
    counts.index = idx.values

    print('Store pancreas atlas to file')
    loompy.create(
        '../data/pancreas_atlas2/dataset.loom',
        layers={'': counts.values},
        row_attrs={'GeneName': counts.index.values},
        col_attrs={'CellID': counts.columns.values, 'cellType': cells['cellType'].values},
        )
