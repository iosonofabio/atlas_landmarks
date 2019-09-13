# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/08/19
content:    Preprocess Smillie et al 2019 (human colon)
'''
import os
import sys
import numpy as np
from scipy.io import mmread
import pandas as pd
import loompy



if __name__ == '__main__':

    fdn = '../data_raw/Smillie_2019/'

    print('Read metadata')
    meta = pd.read_csv(fdn+'all.meta2.txt', sep='\t', index_col=0).iloc[1:]

    print('Read features')
    genes_epi = np.loadtxt(fdn+'Epi.genes.tsv', delimiter='\t', dtype='U50')
    genes_fib = np.loadtxt(fdn+'Fib.genes.tsv', delimiter='\t', dtype='U50')
    genes_imm = np.loadtxt(fdn+'Imm.genes.tsv', delimiter='\t', dtype='U50')
    genes = np.intersect1d(np.intersect1d(genes_epi, genes_fib), genes_imm)

    print('Read cells')
    cells_epi = np.loadtxt(fdn+'Epi.barcodes2.tsv', delimiter='\t', dtype='U50')
    cells_fib = np.loadtxt(fdn+'Fib.barcodes2.tsv', delimiter='\t', dtype='U50')
    cells_imm = np.loadtxt(fdn+'Imm.barcodes2.tsv', delimiter='\t', dtype='U50')

    print('Set output file')
    fdn_out = '../data_full/Smillie_2019/'
    fn_out = fdn_out+'dataset.loom'
    os.makedirs(fdn_out, exist_ok=True)

    print('Read FIB count matrix')
    mat_fib = mmread(fdn+'gene_sorted-Fib.matrix.mtx').astype(np.float32)
    meta_fib = meta.loc[cells_fib]

    print('Exclude sick individuals')
    mat_fib = mat_fib.tocsc()
    ind_c = (meta_fib['Health'] == 'Non-inflamed').values.nonzero()[0]
    mat_fib = mat_fib[:, ind_c]
    meta_fib = meta_fib.iloc[ind_c]

    print('Only take shared genes')
    mat_fib = mat_fib.tocsr()
    ind_g = pd.Index(genes_fib).isin(genes).nonzero()[0]
    mat_fib = mat_fib[ind_g]

    print('Make dense matrix')
    mat_fib = mat_fib.todense()

    print('Bootstrap output loom with with fib')
    with loompy.new(fn_out) as dsl:
        dsl.add_columns(
            layers={'': mat_fib},
            row_attrs={
                'GeneName': genes,
                },
            col_attrs={
                'CellID': meta_fib.index.values,
                'CellType': meta_fib['Cluster'].values,
                'NumberOfGenes': meta_fib['nGene'].astype(int).values,
                'NumberOfUMI': meta_fib['nUMI'].astype(int).values,
                'Subject': meta_fib['Subject'].values,
                'Location': meta_fib['Location'].values,
                'Sample': meta_fib['Sample'].values,
                }
            )
    del mat_fib
    del meta_fib

    print('Read EPI count matrix')
    mat_epi = mmread(fdn+'gene_sorted-Epi.matrix.mtx').astype(np.float32)
    meta_epi = meta.loc[cells_epi]

    print('Exclude sick individuals')
    mat_epi = mat_epi.tocsc()
    ind_c = (meta_epi['Health'] == 'Non-inflamed').values.nonzero()[0]
    mat_epi = mat_epi[:, ind_c]
    meta_epi = meta_epi.iloc[ind_c]

    print('Only take shared genes')
    mat_epi = mat_epi.tocsr()
    ind_g = pd.Index(genes_epi).isin(genes).nonzero()[0]
    mat_epi = mat_epi[ind_g]

    print('Make dense matrix')
    mat_epi = mat_epi.todense()

    print('Add epi to loom file')
    with loompy.connect(fn_out) as dsl:
        dsl.add_columns(
            layers={'': mat_epi},
            col_attrs={
                'CellID': meta_epi.index.values,
                'CellType': meta_epi['Cluster'].values,
                'NumberOfGenes': meta_epi['nGene'].astype(int).values,
                'NumberOfUMI': meta_epi['nUMI'].astype(int).values,
                'Subject': meta_epi['Subject'].values,
                'Location': meta_epi['Location'].values,
                'Sample': meta_epi['Sample'].values,
                }
            )
    del mat_epi
    del meta_epi

    print('Read IMM count matrix')
    mat_imm = mmread(fdn+'gene_sorted-Imm.matrix.mtx').astype(np.float32)
    meta_imm = meta.loc[cells_imm]

    print('Exclude sick individuals')
    mat_imm = mat_imm.tocsc()
    ind_c = (meta_imm['Health'] == 'Non-inflamed').values.nonzero()[0]
    mat_imm = mat_imm[:, ind_c]
    meta_imm = meta_imm.iloc[ind_c]

    print('Only take shared genes')
    mat_imm = mat_imm.tocsr()
    ind_g = pd.Index(genes_imm).isin(genes).nonzero()[0]
    mat_imm = mat_imm[ind_g]

    print('Make dense matrix')
    mat_imm = mat_imm.todense()

    print('Add imm to loom file')
    with loompy.connect(fn_out) as dsl:
        dsl.add_columns(
            layers={'': mat_imm},
            col_attrs={
                'CellID': meta_imm.index.values,
                'cellType': meta_imm['Cluster'].values,
                'NumberOfGenes': meta_imm['nGene'].astype(int).values,
                'NumberOfUMI': meta_imm['nUMI'].astype(int).values,
                'Subject': meta_imm['Subject'].values,
                'Location': meta_imm['Location'].values,
                'Sample': meta_imm['Sample'].values,
                }
            )
    del mat_imm
    del meta_imm
