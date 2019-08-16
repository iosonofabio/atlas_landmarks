# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/07/19
content:    Preprocess pancreas data from Baron et al. 2016.
'''
import os
import sys
import glob
import numpy as np
import pandas as pd
import loompy


if __name__ == '__main__':

    fn_pancreas_atlas_meta = '../data/pancreas_atlas3/GSE84133_family.soft'

    print('Parse brain atlas 3 metadata')
    samples = {'name': [], 'countsUrl': [], 'organism': [], 'disease': []}
    with open(fn_pancreas_atlas_meta, 'rt') as f:
        for line in f:
            if line.startswith('^SAMPLE ='):
                # Pad missing info
                for col in samples:
                    if len(samples[col]) < len(samples['name']):
                        samples[col].append('')
                name = line.rstrip('\n').split()[-1]
                samples['name'].append(name)
            elif line.startswith('!Sample_supplementary_file_1 = '):
                ftp_url = line.rstrip('\n').split()[-1]
                samples['countsUrl'].append(ftp_url)
            elif line.startswith('!Sample_organism_ch1 ='):
                org = line.rstrip('\n').split('=')[-1][1:]
                samples['organism'].append(org)
            elif line.startswith('!Sample_characteristics_ch1 = type 2 diabetes mellitus:'):
                dis = line.rstrip('\n').split(':')[-1][1:]
                samples['disease'].append(dis)
    # Pad missing info
    for col in samples:
        if len(samples[col]) < len(samples['name']):
            samples[col].append('')

    samples = pd.DataFrame(samples).set_index('name')
    samples['countsFilename'] = ['../data/pancreas_atlas3/'+x.split('/')[-1] for x in samples['countsUrl'].values]

    print('Exclude disease')
    ind = samples['disease'] != 'Yes'
    samples = samples.loc[ind]

    print('Keep only human for now')
    ind_human = samples['organism'] == 'Homo sapiens'
    samples_human = samples.loc[ind_human]

    print('Download and parse human pancreas atlas counts')
    from ftplib import FTP
    ncells = len(samples)
    data = []
    for ic, cn in enumerate(samples_human.index):
        print('Loading sample {:}'.format(ic+1))
        fn = samples_human.at[cn, 'countsFilename']
        if not os.path.isfile(fn):
            print('{:}: download file via FTP...'.format(cn), end='')
            url = samples_human.at[cn, 'countsUrl']
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

        datum = pd.read_csv(fn, sep=',', index_col=0)

        meta = datum[['barcode', 'assigned_cluster']]
        counts = datum.iloc[:, 2:-1].T.astype(np.float32)
        ncells = counts.shape[1]
        ngenes = counts.shape[0]
        data.append({
            'meta': meta,
            'counts': counts,
            'ncells': ncells,
            'ngenes': ngenes,
            })

    N = sum(d['ncells'] for d in data)
    L = data[0]['ngenes']
    genes = data[0]['counts'].index
    cells = pd.concat([d['meta'] for d in data])
    matrix = np.empty((L, N), np.float32)
    i = 0
    for d in data:
        n = d['ncells']
        matrix[:, i: i+n] = d['counts'].values
        i += n
    counts = pd.DataFrame(matrix, index=genes, columns=cells.index)
    cells.rename(columns={'assigned_cluster': 'cellType'}, inplace=True)

    print('Store pancreas atlas to file')
    loompy.create(
        '../data/pancreas_atlas3/dataset.loom',
        layers={'': counts.values},
        row_attrs={'GeneName': counts.index.values},
        col_attrs={'CellID': counts.columns.values, 'cellType': cells['cellType'].values},
        )
