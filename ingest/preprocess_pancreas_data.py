# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/07/19
content:    Preprocess pancreas data from Enge et al. 2017.
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    fn_pancreas_atlas_meta = '../data/pancreas_atlas/GSE81547_family.soft'

    print('Parse brain atlas metadata')
    cells = {'name': [], 'cellType': [], 'countsUrl': []}
    with open(fn_pancreas_atlas_meta, 'rt') as f:
        for line in f:
            if line.startswith('^SAMPLE ='):
                name = line.rstrip('\n').split()[-1]
                cells['name'].append(name)
            elif line.startswith('!Sample_characteristics_ch1 = inferred_cell_type:'):
                ctype = line.rstrip('\n').split(':')[-1][1:]
                cells['cellType'].append(ctype)
            elif line.startswith('!Sample_supplementary_file_1 = '):
                ftp_url = line.rstrip('\n').split()[-1]
                cells['countsUrl'].append(ftp_url)
    cells = pd.DataFrame(cells).set_index('name')
    cells['countsFilename'] = ['../data/pancreas_atlas/'+x.split('/')[-1] for x in cells['countsUrl'].values]

    print('Exclude unsure cell types')
    ind = cells['cellType'] != 'unsure'
    cells = cells.loc[ind]

    print('Download and parse pancreas atlas counts')
    from ftplib import FTP
    ncells = len(cells)
    counts = None
    for ic, cn in enumerate(cells.index):
        fn = cells.at[cn, 'countsFilename']
        if not os.path.isfile(fn):
            print('{:}: download file via FTP...'.format(cn), end='')
            url = cells.at[cn, 'countsUrl']
            ftp_root = 'ftp.ncbi.nlm.nih.gov'
            url_fdn = os.path.dirname(url[len('ftp://')+len(ftp_root)+1:])
            url_fn = os.path.basename(url)
            ftp = FTP(ftp_root)
            ftp.login()
            ftp.cwd(url_fdn)
            with open(fn, 'wb') as f:
                ftp.retrbinary(
                    'RETR {:}'.format(url_fn),
                    f.write,
                    )
            try:
                ftp.quit()
            except Exception:
                pass
            print('done!')
        count = pd.read_csv(fn, sep='\t', index_col=0, squeeze=True)
        if counts is None:
            ngenes = count.shape[0]
            genes = count.index.str.rstrip(' ')
            counts = np.zeros((ngenes, ncells), np.float32)
        counts[:, ic] = count.values
    counts = pd.DataFrame(counts, index=genes, columns=cells.index)

    print('Store pancreas atlas to file')
    loompy.create(
        '../data/pancreas_atlas/dataset.loom',
        layers={'': counts.values},
        row_attrs={'GeneName': counts.index.values},
        col_attrs={'CellID': counts.columns.values, 'cellType': cells['cellType'].values},
        )
