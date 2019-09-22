# vim: fdm=indent
'''
author:     Fabio Zanini
date:       15/09/19
content:    Preprocess Young et al 2018 (human kidney)
'''
import os
import sys
import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix
import pandas as pd
import loompy



if __name__ == '__main__':

    fdn = '../data_raw/Young_2018/'

    print('Read metadata')
    meta = pd.read_csv(fdn+'cell_metadata.tsv', sep='\t', index_col='DropletID')

    print('Read features')
    meta_genes = pd.read_csv(fdn+'tableOfCounts_rowLabels.tsv', sep='\t', index_col=1)

    print('Read counts')
    counts = mmread(fdn+'tableOfCounts.mtx').astype(np.float32)

    print('Exclude genes that have two Ensembl IDs')
    tmp = meta_genes['Symbol'].value_counts()
    tmp = tmp.index[tmp > 1]
    ind_uniques = ~(meta_genes['Symbol'].isin(tmp)).values

    meta_genes = meta_genes.loc[ind_uniques]
    counts = counts.tocsr()[ind_uniques].tocsc()

    print('Exclude fetal and tumor')
    ind_healthy = meta['Compartment'].str.startswith('Normal').values
    meta = meta[ind_healthy]
    counts = counts[:, ind_healthy]

    print('Turns out now we can convert to dense')
    counts = counts.todense()

    # FIXME: figure out what those annotations actually mean!!

    print('Set output file')
    fdn_out = '../data_full/Young_2018/'
    fn_out = fdn_out+'dataset.loom'
    os.makedirs(fdn_out, exist_ok=True)

    print('Write to loom')
    with loompy.new(fn_out) as dsl:
        dsl.add_columns(
            layers={'': counts},
            row_attrs={
                'GeneName': meta_genes['Symbol'].values,
                },
            col_attrs={
                'CellID': meta.index.values,
                'CellType': meta['ClusterID'].values,
                'NumberOfGenes': meta['nGenes'].astype(int).values,
                'NumberOfUMI': meta['nUMI'].astype(int).values,
                'Subject': meta['Source'].values,
                'Location': meta['Compartment'].values,
                }
            )
