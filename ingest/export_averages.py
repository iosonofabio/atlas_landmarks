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
import yaml


def discover_datasets():
    '''Discover atlases and list them'''
    fdn_full = '../data_full/'
    datasets = {}
    for sfdn in os.listdir(fdn_full):
        if not os.path.isdir(fdn_full+sfdn):
            continue
        fns = os.listdir(fdn_full+sfdn)
        if len(fns) == 0:
            continue
        dataset = []
        for fn in fns:
            if fn == 'dataset.loom':
                dataset.append(None)
            elif fn.startswith('dataset') and fn.endswith('.loom'):
                tissue = '_'.join(fn.split('.')[0].split('_')[1:])
                dataset.append(tissue)
        if dataset:
            datasets[sfdn] = dataset

    return datasets


class AtlasAverager():
    def __init__(self, name, tissue, overwrite=False):
        self.name = name
        self.tissue = tissue
        self.overwrite = overwrite
        self.get_full_filename()
        self.get_custom_filters()

    def get_full_filename(self):
        fdn = '../data_full/'+self.name
        if self.tissue is None:
            self.full_filename = fdn+'/dataset.loom'
        else:
            self.full_filename = fdn+'/dataset_{:}.loom'.format(self.tissue)

    def get_output_filename(self, metaname):
        fdn = '../data/averages/'
        if self.tissue is None:
            return fdn+metaname+'.loom'
        else:
            return fdn+metaname+'_'+self.tissue+'.loom'

    def get_custom_filters(self):
        '''Some datasets have more than one export'''
        filters = {
                'Darmanis_2015': {'nofetal': lambda x: 'fetal' not in x},
        }
        filters = filters.get(self.name, dict())
        filters[''] = None
        self.filters = filters

    def get_atlas_metadata(self, name=None):
        if name is None:
            name = self.name
        fn = '../atlas_metadata.yml'
        with open(fn, 'rt') as f:
            meta = yaml.safe_load(f)
        return meta.get(name, dict())

    def process_atlas(self):
        print(self.name)

        if self.tissue is None:
            print('Check output files')
        else:
            print('{:}: Check output files'.format(self.tissue))
        fns_out = {}
        for filtname in self.filters:
            if filtname:
                print('Export data, {:}'.format(filtname))
                metaname = self.name+'_'+filtname
            else:
                print('Export data')
                metaname = self.name
            fns_out[filtname] = self.get_output_filename(metaname)

        if (not self.overwrite) and all(os.path.isfile(x) for x in fns_out.values()):
            print('Exists already, skipping')
            return

        if self.tissue is None:
            print('Read data and average by cell type')
        else:
            print('{:}: read data and average by cell type'.format(self.tissue))
        with loompy.connect(self.full_filename) as dsl:
            cts = dsl.ca['cellType']
            n_cells = Counter(cts)
            n_cells_tot = sum(n_cells.values(), 0)
            ctu = list(n_cells.keys())
            N = len(ctu)

            # Exclude ERCC spike-ins and QC features
            features = dsl.ra['GeneName']
            exclude_list = [
                'too_low_aQual',
                'alignment_not_unique',
                'ambiguous',
                'no_feature',
                'not_aligned',
                ]
            ind_fea = []
            for fea in features:
                fea_bool = fea not in exclude_list
                fea_bool &= not fea.startswith('ERCC-')
                fea_bool &= not fea.startswith('_')
                ind_fea.append(fea_bool)
            features = features[ind_fea]
            L = len(features)

            matrix = np.zeros((L, N), dtype=np.float32)
            lstring = max([len(ct) for ct in n_cells])
            cnames = np.array(ctu, dtype='U'+str(lstring + 12))
            ncnames = np.array([n_cells[x] for x in ctu])
            for i, ct in enumerate(ctu):
                print('Cell type: {:}'.format(ct))

                ind_cell = (cts == ct)

                submat = dsl[:, ind_cell]
                submat = submat[ind_fea]
                submat = submat.astype(np.float32)

                # Normalize
                submat *= 1e6 / submat.sum(axis=0)

                # Aritmetic average
                matrix[:, i] = submat.mean(axis=1)

        for filtname, fun in self.filters.items():
            if filtname:
                print('Export data, {:}'.format(filtname))
                metaname = self.name+'_'+filtname
            else:
                print('Export data')
                metaname = self.name

            fn_out = fns_out[filtname]

            file_attrs = self.get_atlas_metadata(metaname)
            file_attrs['NumberOfCells'] = n_cells_tot
            if self.tissue is not None:
                file_attrs['Tissue'] = self.tissue

            # Filter samples and metadata
            if fun is None:
                cnames_filt = cnames
                ncnames_filt = ncnames
                matrix_filt = matrix
            else:
                ind = [fun(x) for x in cnames]
                cnames_filt = cnames[ind]
                ncnames_filt = ncnames[ind]
                matrix_filt = matrix[:, ind]

            loompy.create(
                fn_out,
                layers={'': matrix_filt},
                row_attrs={'GeneName': features},
                col_attrs={
                    'CellType': cnames_filt,
                    'NumberOfCells': ncnames_filt,
                    },
                file_attrs=file_attrs,
                )


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--overwrite', action='store_true')
    args = pa.parse_args()

    datasets = discover_datasets()

    for dsname, tissues in datasets.items():
        for tissue in tissues:
            exporter = AtlasAverager(dsname, tissue, overwrite=args.overwrite)
            exporter.process_atlas()
