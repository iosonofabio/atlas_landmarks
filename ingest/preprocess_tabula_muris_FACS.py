# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/19
content:    Preprocess tabula muris 2018 FACS data.
'''
import os
import sys
import glob
import argparse
from collections import Counter
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    print('Preprocess tabula muris FACS')
    tissues = [
        'Aorta', 'Bladder', 'Brain', 'Diaphragm', 'Fat', 'Heart', 'Kidney',
        'Intestine', 'Muscle', 'Liver', 'Lung', 'Mammary_Gland', 'Marrow',
        'Pancreas', 'Skin', 'Spleen', 'Thymus', 'Tongue', 'Trachea',
        ]
    print('Read metadata')
    fn = '../data_raw/tabula_muris_2018/FACS/tabula_muris_facs_annotations.csv'
    meta = pd.read_csv(fn, sep=',', index_col='cell')['cell_ontology_class']
    for tissue in tissues:
        print(tissue)
        print('Read data and average')
        fns = glob.glob('../data_raw/tabula_muris_2018/FACS/*{:}*.csv.gz'.format(tissue))
        matrices_raw = []
        genes_raw = []
        cells_raw = []
        for ifn, fn in enumerate(fns):
            if len(fns) > 1:
                print('File {:} of {:}'.format(ifn+1, len(fns)))
            m = pd.read_csv(fn, sep=',', compression='gzip', index_col=0)

            # Many don't make it into the metadata file... low quality?
            ind = m.columns.isin(meta.index)
            m = m.loc[:, ind]
            matrices_raw.append(m.values.astype(np.float32))
            genes_raw.append(m.index)
            cells_raw.append(m.columns)

        cells = np.concatenate(cells_raw)
        cts = meta.loc[cells].astype(str)

        # Exclude cells without a cell type
        cts = cts.loc[cts != 'nan']
        ctu = np.unique(cts)
        n_fixed = len(ctu)
        n_cells_d = Counter(cts.values)
        n_cells = np.array([n_cells_d[ct] for ct in ctu])
        n_cells_total = sum(n_cells)

        # Intersection of features
        features = genes_raw[0].values
        if len(genes_raw) > 1:
            for ge in genes_raw[1:]:
                features = np.intersect1d(features, ge.values)
        # Exclude ERCC spike-ins and QC features
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        # Some tissues (e.g. brain) have more than one CSV file for counts
        matrix = np.zeros((L, n_cells_total), dtype=np.float32)
        ni = 0
        for ce, ge, m in zip(cells_raw, genes_raw, matrices_raw):
            # Check that they have a cell type
            ind = ce.isin(cts.index)
            ce = ce[ind]
            m = m[:, ind]

            # Restrict to the right features
            ind_ge = ge.isin(features)
            m = m[ind_ge]

            matrix[:, ni: ni+m.shape[1]] = m

            ni += m.shape[1]

        print('Export full data')
        os.makedirs('../data_full/Tabula_muris_2018_FACS', exist_ok=True)
        loompy.create(
            '../data_full/Tabula_muris_2018_FACS/dataset_{:}.loom'.format(tissue.lower()),
            layers={'': matrix},
            row_attrs={'GeneName': features},
            col_attrs={
                'cellType': cts.values,
                'CellID': cts.index.values,
                },
            )
