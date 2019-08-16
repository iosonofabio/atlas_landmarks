# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/08/19
content:    Export various atlas averages to file for github repo.
'''
import os
import sys
import glob
import argparse
from collections import Counter
import numpy as np
import pandas as pd
import loompy

sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, CountsTable, FeatureSheet, SampleSheet


if __name__ == '__main__':

    print('Baron et al. 2016')
    print('Read data and average')
    with loompy.connect('../data/pancreas_atlas3/dataset.loom') as dsl:
        cts = dsl.ca['cellType']
        ctu = np.unique(cts)
        n_fixed = len(ctu)
        n_cells = Counter(cts)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        # Exclude ERCC spike-ins and QC features
        features = dsl.ra['GeneName']
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ind = cts == ctu[i]
            submat = dsl[:, ind]
            submat = submat[ind_fea]

            # Normalize
            submat *= 1e6 / submat.sum(axis=0)

            # Aritmetic average
            matrix[:, i] = submat.mean(axis=1)

    print('Export data')
    loompy.create(
        '../data/export_averages/human_pancreas_Baron_2016.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'pancreas',
            'Age': 'adult',
            'Number of cells': sum(n_cells),
            },
        )

    sys.exit()

    print('Stitzel et al. 2016')
    print('Read data and average')
    with loompy.connect('../data/pancreas_atlas2/dataset.loom') as dsl:
        cts = dsl.ca['cellType']
        ctu = np.unique(cts)
        n_fixed = len(ctu)
        n_cells = Counter(cts)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        # Exclude ERCC spike-ins and QC features
        features = dsl.ra['GeneName']
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ind = cts == ctu[i]
            submat = dsl[:, ind]
            submat = submat[ind_fea]

            # Normalize
            submat *= 1e6 / submat.sum(axis=0)

            # Aritmetic average
            matrix[:, i] = submat.mean(axis=1)

    print('Export data')
    loompy.create(
        '../data/export_averages/human_pancreas_Stitzel_2016.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'pancreas',
            'Age': 'adult',
            'Number of cells': sum(n_cells),
            },
        )

    print('Croote et al. 2018')
    print('Read data and average')
    with loompy.connect('../data/Bcells_croote/dataset.loom') as dsl:
        cts = dsl.ca['cellType']
        ctu = np.unique(cts)
        n_fixed = len(ctu)
        n_cells = Counter(cts)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        # Exclude ERCC spike-ins and QC features
        features = dsl.ra['GeneName']
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ind = cts == ctu[i]
            submat = dsl[:, ind]
            submat = submat[ind_fea]

            # Normalize
            submat *= 1e6 / submat.sum(axis=0)

            # Aritmetic average
            matrix[:, i] = submat.mean(axis=1)

    print('Export data')
    loompy.create(
        '../data/export_averages/human_Bcells_Croote_2018.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'peripheral blood',
            'Age': 'children',
            'Number of cells': sum(n_cells),
            },
        )

    print('Zanini et al. 2018')
    print('Read data and average')
    with loompy.connect('../data/pbmc_zanini/dataset.loom') as dsl:
        cts = dsl.ca['cellType']
        ctu = list(np.unique(cts))
        ctu.remove('unknown')
        ctu = np.array(ctu)
        n_fixed = len(ctu)
        n_cells = Counter(cts)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        # Exclude ERCC spike-ins and QC features
        features = dsl.ra['GeneName']
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ind = cts == ctu[i]
            submat = dsl[:, ind]
            submat = submat[ind_fea]

            # Data is already normalized

            # Aritmetic average
            matrix[:, i] = submat.mean(axis=1)

    print('Export data')
    loompy.create(
        '../data/export_averages/human_blood_Zanini_2018.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'peripheral blood',
            'Age': 'adult',
            'Number of cells': sum(n_cells),
            },
        )

    print('Enge et al. 2017')
    print('Read data and average')
    with loompy.connect('../data/pancreas_atlas/dataset.loom') as dsl:
        cts = dsl.ca['cellType']
        ctu = np.unique(cts)
        n_fixed = len(ctu)
        n_cells = Counter(cts)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        # Exclude ERCC spike-ins and QC features
        features = dsl.ra['GeneName']
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ind = cts == ctu[i]
            submat = dsl[:, ind]
            submat = submat[ind_fea]

            # Normalize (for averaging)
            submat *= 1e6 / submat.sum(axis=0)

            # Aritmetic average with pseudocounts
            matrix[:, i] = submat.mean(axis=1)

    print('Export data')
    loompy.create(
        '../data/export_averages/human_pancreas_Enge_2017.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'pancreas',
            'Age': 'fetal and adult',
            'Number of cells': sum(n_cells),
            },
        )

    print('Tabula muris FACS')
    tissues = [
        'Aorta', 'Bladder', 'Brain', 'Diaphragm', 'Fat', 'Heart', 'Kidney',
        'Intestine', 'Muscle', 'Liver', 'Lung', 'Mammary_Gland', 'Marrow',
        'Pancreas', 'Skin', 'Spleen', 'Thymus', 'Tongue', 'Trachea',
        ]
    print('Read metadata')
    fn = '../data/tabula_muris/tabula_muris_facs_annotations.csv'
    meta = pd.read_csv(fn, sep=',', index_col='cell')['cell_ontology_class']
    for tissue in tissues:
        print(tissue)
        print('Read data and average')
        fns = glob.glob('../data/tabula_muris/FACS/*{:}*.csv.gz'.format(tissue))
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
        n_cells = Counter(cts.values)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        # Intersection of features
        features = genes_raw[0].values
        if len(genes_raw) > 1:
            for ge in genes_raw[1:]:
                features = np.intersect1d(features, ge.values)
        # Exclude ERCC spike-ins and QC features
        ind_fea = [not (fea.startswith('ERCC-') or fea.startswith('_')) for fea in features]
        features = features[ind_fea]
        L = len(features)

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ni = 0
            for ce, ge, m in zip(cells_raw, genes_raw, matrices_raw):
                # Check that they have a cell type
                ind = ce.isin(cts.index)
                ce = ce[ind]
                m = m[:, ind]

                # Restrict to the right cell type
                ctsi = cts.loc[ce]
                ind = ctsi == ctu[i]
                m = m[:, ind]

                # Restrict to the right features
                ind_ge = ge.isin(features)
                submat = m[ind_ge]

                # Normalize (for averaging)
                submat *= 1e6 / submat.sum(axis=0)

                # Aritmetic average
                matrix[:, i] += submat.sum(axis=1)
                ni += ind.sum()

            matrix[:, i] /= ni

        print('Export data')
        loompy.create(
            '../data/export_averages/mouse_{:}_TabulaMuris_2018_FACS.loom'.format(tissue.lower()),
            layers={'': matrix},
            row_attrs={'GeneName': features},
            col_attrs={
                'CellType': ctu,
                'NumberOfCells': n_cells,
                },
            file_attrs={
                'Normalization': 'counts per million reads',
                'Averaging': 'aritmetic after normalization',
                'Species': 'Mus Musculus',
                'Tissue': tissue.lower(),
                'Age': '3 months',
                'Number of cells': sum(n_cells),
                },
            )

    print('Darmanis et al. 2015')
    print('Read data and average')
    with loompy.connect('../data/brain_atlas/dataset.loom') as dsl:
        cts = dsl.ca['cellType']
        features = dsl.ra['GeneName']
        ctu = np.unique(cts)
        n_fixed = len(ctu)
        L = len(features)
        n_cells = Counter(cts)
        n_cells = np.array([n_cells[ct] for ct in ctu])

        matrix = np.zeros((L, n_fixed), dtype=np.float32)
        for i in range(n_fixed):
            ind = cts == ctu[i]
            submat = dsl[:, ind]

            # Normalize (for averaging)
            submat *= 1e6 / submat.sum(axis=0)

            # Aritmetic average with pseudocounts
            matrix[:, i] = submat.mean(axis=1)

    print('Export data (with fetal)')
    loompy.create(
        '../data/export_averages/human_brain_Darmanis_2015.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'brain',
            'Age': 'fetal and adult',
            'Number of cells': sum(n_cells),
            },
        )

    print('Export data (without fetal)')
    ind = np.array(['fetal' not in x for x in ctu])
    matrix = matrix[:, ind]
    ctu = ctu[ind]
    n_cells = n_cells[ind]
    loompy.create(
        '../data/export_averages/human_brain_Darmanis_2015_nofetal.loom',
        layers={'': matrix},
        row_attrs={'GeneName': features},
        col_attrs={
            'CellType': ctu,
            'NumberOfCells': n_cells,
            },
        file_attrs={
            'Normalization': 'counts per million reads',
            'Averaging': 'aritmetic after normalization',
            'Species': 'Homo Sapiens',
            'Tissue': 'brain',
            'Age': 'adult only (no fetal)',
            'Number of cells': sum(n_cells),
            },
        )
