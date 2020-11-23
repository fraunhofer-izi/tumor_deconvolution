#! /usr/bin/env python
#SBATCH --job-name="generate cell type characteristics"
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=10-00:00:00
#SBATCH --output=logs/genchar.%j.slurmlog

import warnings
import os
import sys
import numpy as np
if __name__ == "__main__" and not sys.argv[-1].endswith('json'):
    import argparse
    desc = 'Makes cell type characteristics.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('grain', help='"coarse" or "fine" grain characteristics',
                        metavar='grain', type=str, nargs='?', default='coarse')
    parser.add_argument('inFolder', help='Folder with expressions in tsv files per grain.',
                        metavar='inFolder',type=str, nargs='?', default='../data/tissue_tsv_HUGO_array/')
    parser.add_argument('outFile', help='File in which the characteristics saved.',
                        metavar='outFile',type=str, nargs='?', default='../data/genchar_coarse_unclean_array.pkl')
    parser.add_argument('--nSampls', help='Number of Dirichlet samples for incremental PCA (default is 1e6).',
            type=int, default=int(1e6), metavar='integer')
    parser.add_argument('--center', help='Center whitening per cell type.',
            action='store_true')
    parser.add_argument('--simp', help='Do second pca only on cell type means (effect only with --centered).',
            action='store_true')
    parser.add_argument('--nComp', help='Number of components for non-simple characterization (default is 100).',
            type=int, default=100, metavar='integer')
    parser.add_argument('--dcnComp', help='Number of components for incrementel PCA on dirichlet samples (default is 1000).',
            type=int, default=1000, metavar='integer')
    parser.add_argument('--minquality', help='Minimum value of aggregated_qaulity for sample to be included.',
            type=np.float, default=6, metavar='float')
    parser.add_argument('--nonanquality', help='Do not include samples without a quality tag.',
            action='store_true')
    args = parser.parse_args()
    warnings.filterwarnings('ignore', category=FutureWarning)
elif sys.argv[-1].endswith('json'):
    # We are probably in a jupyter notebook
    class A:
        center = True
        simp = True
        grain = 'fine'
        nonanquality = False
        nComp = 100
        minquality = 6
        nSampls = 2e4
        inFolder = '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/tissue_tsv_whseq_array/'
        outFile = '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/genchar_coarse_unclean_center_simp_whseq_interactive.pkl'
    args = A()
    if os.getcwd().endswith('python'):
        os.chdir('..')
import pickle
import re
import tempfile
from collections import Counter
from itertools import tee, count
from tqdm.auto import tqdm
import random
from scipy.stats import dirichlet, multivariate_normal
from scipy.special import digamma
from sklearn.decomposition import IncrementalPCA, PCA
import plotnine as pn
import pandas as pd
import multiprocess as mp
import tempfile
from property_deps import property_deps
import pymanopt as pman
from pymanopt.solvers import SteepestDescent
from pymanopt.manifolds import Stiefel
import tensorflow as tf


class SuppressPrints:
    def __init__(self, suppress=True):
        self.suppress = suppress

    def __enter__(self):
        if self.suppress:
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.suppress:
            sys.stdout.close()
            sys.stdout = self._original_stdout


def tau_inv(log_draw):
    o_samp = log_draw - (log_draw.sum(axis=1, keepdims=True) / log_draw.shape[1])
    return o_samp

def uniquify(seq, suffs = count(1)):
    not_unique = [k for k,v in Counter(seq).items() if v>1]
    suff_gens = dict(zip(not_unique, tee(suffs, len(not_unique))))
    for s in seq:
        try:
            suffix = str(next(suff_gens[s]))
        except KeyError:
            yield s
            continue
        else:
            yield s + '.' + suffix

class Expressions:

    def __init__(self, grain='coarse', TSV_dir='../data/tissue_tsv/',
                 anno_file='../data/sample_data.tsv', exclude_test_samples=False):
        self.grain = grain
        self.TSV_dir = TSV_dir
        self.anno_file = anno_file
        self.exclude_test_samples = exclude_test_samples

    @property_deps('alphas', 'test_counts')
    def expression_tsv(self):
        expression_tsv = dict()
        grain = self.grain
        TSV_dir = self.TSV_dir
        for fl in os.listdir(TSV_dir):
            if re.match(r'{}-.*.tsv'.format(grain), fl) is not None \
            and re.match(r'.*too.unspecific.tsv', fl) is None:
                cell_type = re.sub(r'{}-(.*).tsv'.format(grain), r'\1', fl)
                expression_tsv[cell_type] = TSV_dir + fl
        return expression_tsv

    @property_deps()
    def cell_types(self):
        return list(self.alphas.keys())

    @property_deps('test_counts')
    def anno(self):
        return pd.read_csv(self.anno_file, sep = "\t", low_memory=False).set_index('id')

    def _set_alphas(self):
        manager = mp.Manager()
        shared_alphas = manager.dict()
        shared_test_counts = manager.dict()
        self.alphas = dict()
        self.test_counts = dict()
        test_anno = self.anno.groupby('id')['test_data'].agg(is_test=any)
        exclude_test_samples = self.exclude_test_samples
        def load_alphas(entry):
            data = pd.read_csv(entry[1], sep = "\t", skip_blank_lines=False,
                               keep_default_na=False)
            samp_names = [re.sub('\.[0-9]+', '', x) for x in data.columns]
            if not samp_names:
                return entry[0], False
            is_test_sample = [test_anno.loc[x, 'is_test'] if x in test_anno.index
                              else False for x in samp_names]
            is_test_sample = np.array(is_test_sample)
            if exclude_test_samples is True and any(is_test_sample):
                shared_alphas[entry[0]] = data.loc[:,~is_test_sample] + 1
                shared_test_counts[entry[0]] = data.loc[:,is_test_sample]
            else:
                try:
                    shared_alphas[entry[0]] = data + 1
                except Exception as e:
                    raise Exception(f'Unable to deal with {entry[1]}')
            return entry[0], True
        tasks = list(self.expression_tsv.items())
        with mp.Pool() as pool:
            for key, success in tqdm(pool.imap(load_alphas, tasks),
                            total=len(tasks), desc='loading counts'):
                if success is False:
                    warnings.warn('No samples found for {}.'.format(key), RuntimeWarning)
                    continue
                self.alphas[key] = shared_alphas[key]
                if key in shared_test_counts.keys():
                    self.test_counts[key] = shared_test_counts[key]
        return self.alphas, self.test_counts

    @property_deps('log_E_p', 'E_log_p', 'features', 'cell_types')
    def alphas(self):
        alphas, _ = self._set_alphas()
        return alphas

    @property_deps()
    def test_counts(self):
        _, test_counts = self._set_alphas()
        return test_counts

    @property_deps()
    def features(self):
        features = self.alphas[self.cell_types[0]].index
        return list(features)

    @property_deps('means', 'tau_projected_values')
    def log_E_p(self):
        _ = self.alphas # trigger lazy eval for better tqdm
        def get_expected_p(key):
            l_alphas = self.alphas[key].values.T
            return np.log(l_alphas) - np.log(np.sum(l_alphas, axis=1, keepdims=True))
        return [get_expected_p(ct) for ct in tqdm(self.cell_types, 'log_E_p')]

    @property_deps()
    def E_log_p(self):
        def get_expected_p(key):
            l_alphas = self.alphas[key].values.T
            return digamma(l_alphas) - digamma(np.sum(l_alphas, axis=1, keepdims=True))
        return [get_expected_p(ct) for ct in tqdm(self.cell_types, 'E_log_p')]

    @property_deps('named_means')
    def means(self):
        return [np.mean(x, axis=0) for x in self.tau_projected_values]

    @property_deps()
    def named_means(self):
        return dict(zip(self.cell_types, self.means))

    @property_deps('transformed_expr')
    def tau_projected_values(self):
        pvals = list()
        for logep in tqdm(self.log_E_p, desc='tau projection'):
            pvals.append(tau_inv(logep))
        return pvals

    def get_univariate_dist(self):
        values = np.concatenate(self.tau_projected_values, axis=0)
        gene_means = np.mean(values, axis=0)
        gene_sds = np.std(values, axis=0)
        return gene_means, gene_sds

class ResampledExperssions(Expressions):

    def __init__(self, expr, n_per_type=1e4, seed=7, tmp_dir='/dev/shm/',
                 blacklist=list(), ncores=1):
        super(ResampledExperssions, self).__init__(grain=expr.grain,
                                                   TSV_dir=expr.TSV_dir,
                                                   anno_file=expr.anno_file)
        self.expr = expr
        self.test_counts = expr.test_counts
        self.n_per_type = n_per_type
        self.seed = seed
        self.tmp_dir = tmp_dir
        self.blacklist = blacklist
        self.ncores = ncores
        self.anno = expr.anno

    def _set_alphas(self):
        shared_alphas = self.manager.dict()
        in_alphas = self.selected_alphas
        self.alphas = dict()
        n = int(self.n_per_type)
        np.random.seed(seed=self.seed)
        def draw_samples(key):
            df = in_alphas[key]
            seed = hash(key) % 2**32
            np.random.seed(seed=seed)
            samples = df.T.sample(n=n, random_state=seed+1, replace=True)
            names = list()
            draws = list()
            for name, samp in samples.iterrows():
                draws.append(dirichlet.rvs(samp, size=1).flatten())
                names.append(name)
            names = uniquify(names)
            new_df = pd.DataFrame(dict(zip(names, draws)), index=df.index)
            shared_alphas[key] = new_df
            return key
        if self.ncores == 1:
            for key in in_alphas.keys():
                print('Doing {} ...'.format(key))
                _ = draw_samples(key)
                self.alphas[key] = shared_alphas[key]
        else:
            with mp.Pool(self.ncores) as pool:
                for key in tqdm(pool.imap(draw_samples, in_alphas.keys()),
                                total=len(in_alphas), desc='generating counts'):
                    self.alphas[key] = shared_alphas[key]
        return self.alphas, None

    def resample(self):
        del self.alphas

    @property_deps('alphas', 'selected_alphas')
    def manager(self):
        return mp.Manager()

    @property_deps('means', 'tau_projected_values')
    def log_E_p(self):
        return [np.log(v.T.values) for v in self.alphas.values()]

    @property
    def E_log_p(self):
        return self.log_E_p

    @property_deps('alphas')
    def selected_alphas(self):
        result = self.manager.dict()
        for key, values in self.expr.alphas.items():
            if key in self.background:
                continue
            ind = [samp not in self.blacklist for samp in values.columns]
            if any(ind):
                result[key] = values.loc[:, ind]
        return result

    @property_deps('selected_alphas')
    def background(self):
        return set(['CRC', 'crc', 'BRCA', 'brca', 'others'])

class Dirichlet_PCA:

    def __init__(self, expr, samples=1e6, n_components=1000, seed=7):
        assert isinstance(expr, Expressions), 'experssions musst be of class Expressions'
        self.expr = expr
        self.total_samples = int(samples)
        if n_components is not None:
            self.n_components = n_components
        self.seed = seed

    @property_deps('transformer')
    def n_components(self):
        return np.min(self.expr.shape)

    @property_deps('transformed_expr', 'is_fitted')
    def transformer(self):
        return IncrementalPCA(whiten=True, copy=False,
                              n_components=self.n_components)

    def make_sample(self, n, seed=None, pool=None):
        if seed is None:
            seed = self.seed
        if pool is None:
            pool = mp.Pool()
            close_pool = True
        else:
            close_pool = False
        def make_dirichlet(task):
            seed, sample = task
            draw = dirichlet.rvs(sample.T[0], size=1, random_state=seed+1)
            return tau_inv(np.log(draw))
        random.seed(seed)
        np.random.seed(seed=seed)
        seeds = np.random.randint(1e6, size=n)
        ex = self.expr.alphas
        cts = self.expr.cell_types
        samples = [ex[ct].sample(n=1, axis=1).values for ct in random.sample(cts, 1) for i in range(n)]
        tasks = zip(seeds, samples)
        samples = list(pool.imap_unordered(make_dirichlet, tasks, chunksize=10))
        if close_pool is True:
            pool.close()
        return np.concatenate(samples)

    def fit_transformer(self, iterations=None, samps_per_iteration=None):
        if samps_per_iteration is None:
            samps_per_iteration = 10*self.n_components
        if iterations is None:
            iterations = np.ceil(self.total_samples/samps_per_iteration).astype(int)
        _ = self.expr.alphas # trigger lazy eval
        self.pool = mp.Pool()
        for i in tqdm(range(iterations), desc='incremental pca'):
            samples = self.make_sample(samps_per_iteration,
                                       seed=i+self.seed, pool=self.pool)
            self.transformer.partial_fit(samples)
        del self.transformed_expr
        self.pool.close()
        del self.pool
        self.is_fitted = True
        return self.transformer

    @property_deps('dataframe')
    def transformed_expr(self):
        pvals = np.concatenate(self.expr.tau_projected_values, axis=0)
        return self.transform(pvals)

    @property_deps()
    def dataframe(self):
        df = pd.DataFrame(self.transformed_expr)
        df.columns = ['PC ' + str(i+1) for i in range(self.n_components)]
        cell_type = list()
        for key, value in self.expr.alphas.items():
            cell_type = cell_type + (value.shape[1]*[key])
        df.loc[:, 'cell type'] = cell_type
        return df

    def plot(self, pc1=1, pc2=2, selection=None, color='cell type',
             filter_on='cell type', alpha=.8):
        data = self.dataframe
        if selection is not None:
            ind = [x in selection for x in data[filter_on]]
            data = data.loc[ind,]
        pl = pn.ggplot(pn.aes('PC '+str(pc1), 'PC '+str(pc2), color=color), data) + pn.geom_point(alpha=alpha)
        return pl

    def transform(self, values):
        if not hasattr(self, 'is_fitted') or self.is_fitted is not True:
            self.fit_transformer()
        return self.transformer.transform(values)


class Dirichlet_PCA_centered(Dirichlet_PCA):

    def make_sample(self, n, centers, seed=None, pool=None):
        if seed is None:
            seed = self.seed
        if pool is None:
            pool = mp.Pool()
            close_pool = True
        else:
            close_pool = False
        def make_dirichlet(task):
            seed, sample, center = task
            draw = dirichlet.rvs(sample.T[0], size=1, random_state=seed+1)
            return tau_inv(np.log(draw)) - center
        random.seed(seed)
        np.random.seed(seed=seed)
        seeds = np.random.randint(1e6, size=n)
        ex = self.expr.alphas
        cts = random.choices(self.expr.cell_types, k=n)
        samples = [ex[ct].sample(n=1, axis=1).values for ct in cts]
        cents = [centers[ct] for ct in cts]
        tasks = zip(seeds, samples, cents)
        samples = list(pool.imap_unordered(make_dirichlet, tasks, chunksize=10))
        if close_pool is True:
            pool.close()
        return np.concatenate(samples)

    def fit_transformer(self, iterations=None, samps_per_iteration=None):
        if samps_per_iteration is None:
            samps_per_iteration = 10*self.n_components
        if iterations is None:
            iterations = np.ceil(self.total_samples/samps_per_iteration).astype(int)
        cc = np.mean(np.stack(self.expr.means), axis=0)
        centers = {key:values-cc for key, values in self.expr.named_means.items()}
        self.pool = mp.Pool()
        for i in tqdm(range(iterations), desc='incremental pca'):
            samples = self.make_sample(samps_per_iteration, centers,
                                       seed=i+self.seed, pool=self.pool)
            self.transformer.partial_fit(samples)
        del self.transformed_expr
        self.pool.close()
        del self.pool
        self.is_fitted = True
        return self.transformer


class MultiCoCharacterizer:

    def __init__(self, expressions, n_components=100, seed=8,
                 n_batches=100, batch_size_factor=10):
        assert isinstance(expressions, Expressions), 'experssions musst be of class Expressions'
        self.expr = expressions
        self.n_components = n_components
        self.seed = seed
        self.n_batches = n_batches
        self.batch_size_factor = batch_size_factor

    @property_deps('affine_maps')
    def expr(self):
        return None

    @property_deps('selection')
    def background(self):
        return set(['CRC', 'crc', 'BRCA', 'brca', 'others'])

    @property_deps()
    def background_chars(self):
        bgc = dict()
        expr = np.concatenate(self.expr.tau_projected_values, axis=0)
        bgc['all'] = {'mean':np.mean(expr, axis=0),
                     'std':np.std(expr, axis=0)}
        ptvalues = dict(zip(self.expr.cell_types, self.expr.tau_projected_values))
        for key in self.background:
            if key not in ptvalues:
                continue
            expr = ptvalues[key]
            bgc[key] = {'mean':np.mean(expr, axis=0),
                       'std':np.std(expr, axis=0)}
        return bgc

    @property_deps('affine_maps')
    def selection(self):
        return set(self.expr.cell_types) - self.background

    @property_deps('affine_maps')
    def n_components(self):
        return None

    @property
    def blacklisted(self):
        if 'blacklisted' in self.data.columns:
            return self.data.loc[:, 'blacklisted']
        return pd.Series([False]*self.data.shape[0], index=self.data.index)

    @blacklisted.setter
    def blacklisted(self, values):
        if any(self.blacklisted != values):
            del self.pre_means
            self.data['blacklisted'] = values

    @property_deps()
    def data(self):
        cell_type = list()
        ids = list()
        alphas = self.expr.alphas
        for key, values in alphas.items():
            cell_type += list(values.shape[1]*[key])
            ids += list(values.columns.values)
        data = pd.DataFrame({'id':ids, 'cell type':cell_type}).set_index('id')
        data['samp_name'] = [re.sub('\.[0-9]+', '', x) for x in data['id']]
        anno = self.expr.anno.copy().reset_index()
        anno['samp_name'] = [re.sub('\.[0-9]+', '', x) for x in anno['id']]
        del anno['id']
        data = pd.merge(data, anno, how='left', left_on=['samp_name', 'cell type'],
                        right_on=['samp_name', ch.expr.grain+'.cell.type'])
        return data

    @property
    def batch_size(self):
        return self.n_components * self.batch_size_factor

    @property_deps('affine_trans')
    def pcas(self):
        alphas = self.expr.alphas
        blacklisted = self.blacklisted
        print('Calculating coordinate systems...')
        if self.seed is not None:
            np.random.seed(seed=self.seed)
        def get_transf(cell_type):
            sample_names = alphas[cell_type].columns
            ind = ~blacklisted[sample_names]
            values = alphas[cell_type].loc[:, ind]
            n = sum(ind)
            pca = IncrementalPCA(whiten=True, copy=False, n_components=self.n_components)
            for i in tqdm(range(self.n_batches), desc=cell_type):
                subsample = np.random.randint(n, size=self.batch_size)
                samples = [dirichlet.rvs(values.iloc[:, i], size=1) for i in subsample]
                samples = tau_inv(np.log(np.concatenate(samples)))
                pca.partial_fit(samples)
            return cell_type, pca
        transformations = dict()
        for sel in self.selection:
            key, pca = get_transf(sel)
            transformations[key] = pca
        return transformations

    @property_deps('remaining_variation')
    def affine_trans(self):
        transformations = dict()
        for key, pca in self.pcas.items():
            b = pca.mean_
            Ainv = pca.components_
            A = Ainv.T
            Ainv = np.dot(np.diag(np.sqrt(pca.explained_variance_)), Ainv)
            A = np.dot(A, np.diag(1/np.sqrt(pca.explained_variance_)))
            transformations[key] = dict({'A':A, 'Ainv':Ainv, 'b':b})
        return transformations

    @staticmethod
    def get_variation(values, A, Ainv, b):
        schifted_values = values - b
        co_projection = np.dot(np.dot(schifted_values, A), Ainv)
        schifted_values -= co_projection
        mean = np.mean(schifted_values, axis=0)
        std = np.std(schifted_values, axis=0)
        return mean, std

    @property_deps()
    def remaining_variation(self):
        cell_types = self.expr.cell_types
        blacklisted = self.blacklisted
        ptvalues = dict(zip(cell_types, self.expr.tau_projected_values))
        sample_names = {ct:list(a.columns) for ct, a in self.expr.alphas.items()}
        transformations = self.affine_trans
        def get_variation(cell_type):
            tf = transformations[cell_type]
            A = tf['A']
            Ainv = tf['Ainv']
            b = tf['b']
            ind = ~blacklisted[sample_names[cell_type]]
            values = ptvalues[cell_type][ind, :]
            mean, std = self.get_variation(values, A, Ainv, b)
            return cell_type, mean, std
        variations = dict()
        for sel in tqdm(self.selection, desc='variations'):
            key, mean, std = get_variation(sel)
            variations[key] = dict({'mean':mean, 'std':std})
        return variations


class Characterizer:

    def __init__(self, dirichlet_pca, expressions=None):
        assert isinstance(dirichlet_pca, Dirichlet_PCA), 'dirichlet_pca musst be of class Dirichlet_PCA'
        self.dpca = dirichlet_pca
        if expressions is None:
            self.expr = dirichlet_pca.expr
            self._new_expr = False
        else:
            assert isinstance(expressions, Expressions), 'experssions musst be of class Expressions'
            self.expr = expressions
            self._new_expr = True

    @property_deps('selection')
    def background(self):
        return set(['CRC', 'crc', 'BRCA', 'brca', 'others'])

    @property_deps()
    def background_chars(self):
        bgc = dict()
        expr = np.concatenate(self.expr.tau_projected_values, axis=0)
        bgc['all'] = {'mean':np.mean(expr, axis=0),
                     'std':np.std(expr, axis=0)}
        ptvalues = dict(zip(self.expr.cell_types, self.expr.tau_projected_values))
        for key in self.background:
            if key not in ptvalues:
                continue
            expr = ptvalues[key]
            bgc[key] = {'mean':np.mean(expr, axis=0),
                       'std':np.std(expr, axis=0)}
        return bgc

    @property_deps('transformed_expr', 'is_fitted', 'affine_trans')
    def transformer(self):
        return PCA(copy=True, n_components=self.pre_means.shape[0]-1)

    def fit_transformer(self):
        trans = self.transformer.fit(self.pre_means)
        self.is_fitted = True
        return trans

    @property
    def blacklisted(self):
        if 'blacklisted' in self.data.columns:
            return self.data.loc[:, 'blacklisted']
        return pd.Series([False]*self.data.shape[0], index=self.data.index)

    @blacklisted.setter
    def blacklisted(self, values):
        if any(self.blacklisted != values):
            del self.pre_means
            self.data['blacklisted'] = values
            cell_types = self.expr.cell_types
            for ct in self.expr.cell_types:
                if ct not in self.selection:
                    continue
                samples = self.expr.alphas[ct].columns
                is_blacklisted = values[samples]
                if is_blacklisted.all():
                    raise Exception(f'All {ct} are blacklisted.')

    @property
    def blacklist(self):
        return list(self.blacklisted.index[self.blacklisted])

    @property_deps('pre_means', 'chars')
    def selection(self):
        return set(self.expr.cell_types) - self.background

    @property_deps('pre_means', 'affine_trans')
    def dpca(self):
        return None

    @property_deps('selection', 'data', 'pre_means', 'background_chars')
    def expr(self):
        return None

    @property_deps('transformer')
    def pre_means(self):
        means = list()
        proj_values = self.expr.tau_projected_values
        cell_types = self.expr.cell_types
        for values, ct in zip(proj_values, cell_types):
            if ct not in self.selection:
                continue
            samples = self.expr.alphas[ct].columns
            is_blacklisted = self.blacklisted[samples]
            l_values = values[~is_blacklisted, :]
            means.append(np.mean(l_values, axis=0))
        means = self.dpca.transform(np.stack(means))
        return means

    @property_deps('coords', 'chars')
    def transformed_expr(self):
        if self._new_expr is True:
            pvals = np.concatenate(self.expr.tau_projected_values, axis=0)
            dpcat = self.dpca.transform(pvals)
        else:
            dpcat = self.dpca.transformed_expr
        return self.transform(dpcat)

    @property_deps('coords')
    def data(self):
        cell_type = list()
        ids = list()
        alphas = self.expr.alphas
        for key, values in alphas.items():
            cell_type += list(values.shape[1]*[key])
            ids += list(values.columns.values)
        data = pd.DataFrame({'id':ids, 'cell type':cell_type}).set_index('id')
        anno = self.expr.anno
        re_anno = anno.loc[~anno.index.duplicated(), :].reindex(data.index)
        data['samp_name'] = [re.sub('\.[0-9]+', '', x) for x in data.index]
        data = pd.merge(data, re_anno, how='left', left_on='samp_name', right_index=True)
        return data

    @property_deps()
    def coords(self):
        pDat = pd.DataFrame(self.transformed_expr)
        pDat.columns = ['PC ' + str(i+1) for i in range(self.transformed_expr.shape[1])]
        pDat['id'] = self.data.index
        return pDat.set_index('id')

    @property
    def plot_data(self):
        return pd.concat([self.data, self.coords], axis=1)

    def plot(self, pc1=1, pc2=2, selection=None, color='cell type', shape=None,
             filter_on='cell type', alpha=.8, bl_rm=False, data=None):
        if data is None:
            data = self.plot_data
        if selection is not None:
            ind = [x in selection for x in data[filter_on]]
            data = data.loc[ind,]
        if bl_rm is True:
            data = data.loc[~self.blacklisted,]
        pl = pn.ggplot(pn.aes('PC '+str(pc1), 'PC '+str(pc2), color=color), data) + pn.geom_point(alpha=alpha)
        if shape is not None:
            pl = pl + pn.aes(shape=shape)
        return pl

    @property_deps()
    def chars(self):
        data = self.data
        transformed_expr = self.transformed_expr
        blacklisted = self.blacklisted
        def fit_normal(value, field='cell type'):
            ind = (data[field] == value) & (~blacklisted)
            dat = transformed_expr[ind, :]
            mean = np.mean(dat, axis=0)
            sigma = np.cov(dat, rowvar=False)
            return value, mean, sigma
        sel = self.selection
        chars = dict()
        with mp.Pool(4) as pool:
            for key, mean, sigma in tqdm(pool.imap_unordered(fit_normal, sel),
                                         total=len(sel), desc='character'):
                chars[key] = dict({'mean':mean, 'sigma':sigma})
        return chars

    def calculate_probs(self):
        coords = self.transformed_expr
        ct_map = dict()
        for key in self.chars.keys():
            mean = self.chars[key]['mean']
            sigma = self.chars[key]['sigma']
            try:
                self.data[key + ' prob'] = multivariate_normal.pdf(coords, mean=mean, cov=sigma)
                self.data[key + ' log-prob'] = multivariate_normal.logpdf(coords, mean=mean, cov=sigma)
            except np.linalg.LinAlgError:
                self.data[key + ' prob'] = np.nan
                self.data[key + ' log-prob'] = np.nan
            ct_map[key + ' log-prob'] = key
        self.data['max prop'] = self.data[ct_map.keys()].idxmax(axis=1).map(ct_map)
        self.data['inconsistent'] = self.data['max prop'] != self.data['cell type']
        return

    def transform(self, values):
        if not hasattr(self, 'is_fitted') or self.is_fitted is not True:
            self.fit_transformer()
        return self.transformer.transform(values)

    @property_deps()
    def affine_trans(self):
        if not hasattr(self.dpca, 'is_fitted') or self.dpca.is_fitted is not True:
            self.dpca.fit_transformer()
        t1 = self.dpca.transformer
        if not hasattr(self, 'is_fitted') or self.is_fitted is not True:
            self.fit_transformer()
        t2 = self.transformer
        b1 = t1.mean_
        A1inv = t1.components_
        A1 = A1inv.T
        if t1.whiten:
            A1 = np.dot(A1, np.diag(1/np.sqrt(t1.explained_variance_)))
            A1inv = np.dot(np.diag(np.sqrt(t1.explained_variance_)), A1inv)
        b2 = t2.mean_
        A2inv = t2.components_
        A2 = A2inv.T
        if t2.whiten:
            A2 = np.dot(A2, np.diag(1/np.sqrt(t2.explained_variance_)))
            A2inv = np.dot(np.diag(np.sqrt(t2.explained_variance_)), A2inv)
        b = b1 + np.dot(b2, A1inv)
        A = np.dot(A1, A2)
        Ainv = np.dot(A2inv, A1inv)
        return b, A, Ainv

    def get_remaining_variation(self):
        cell_types = self.data['cell type']
        blacklisted = self.blacklisted
        b, A, Ainv = self.affine_trans
        values = np.concatenate(self.expr.tau_projected_values, axis=0) - b
        co_projection = np.dot(np.dot(values, A), Ainv)
        values -= co_projection
        def get_variation(cell_type):
            ind = (cell_types == cell_type) & (~blacklisted)
            dat = values[ind, :]
            mean = np.mean(dat, axis=0)
            std = np.std(dat, axis=0)
            return cell_type, mean, std
        variations = dict()
        for sel in tqdm(self.selection, desc='variatios'):
            key, mean, std = get_variation(sel)
            variations[key] = dict({'mean':mean, 'std':std})
        return variations

    def get_confusion_matrix(self):
        if 'max prop' not in self.data.keys():
            self.calculate_probs()
        mat = pd.crosstab(self.data['max prop'], self.data['cell type'],
                          rownames=['assigned'], colnames=['label'])
        return mat

class Characterizer_centered(Characterizer):

    def __init__(self, dirichlet_pca_centered, expressions=None, samples=1e6,
                 n_components=100, seed=42):
        msg = 'dirichlet_pca_centered musst be of class Dirichlet_PCA_centere'
        assert isinstance(dirichlet_pca_centered, Dirichlet_PCA_centered), msg
        super().__init__(dirichlet_pca_centered, expressions)
        self.n_components = n_components
        self.total_samples = int(samples)
        self.seed = seed

    @property_deps('transformer')
    def n_components(self):
        return np.min(self.dpca.transformed_expr.shape)

    @property_deps('transformed_expr', 'is_fitted', 'affine_trans')
    def transformer(self):
        return IncrementalPCA(whiten=True, copy=False,
                              n_components=self.n_components)

    def make_sample(self, n, seed=None, pool=None):
        if seed is None:
            seed = self.seed
        if pool is None:
            pool = mp.Pool()
            close_pool = True
        else:
            close_pool = False
        def make_dirichlet(task):
            seed, sample = task
            draw = dirichlet.rvs(sample.T[0], size=1, random_state=seed+1)
            return tau_inv(np.log(draw))
        random.seed(seed)
        np.random.seed(seed=seed)
        seeds = np.random.randint(1e6, size=n)
        ex = self.expr.alphas
        cts = random.choices(self.expr.cell_types, k=n)
        samples = [ex[ct].sample(n=1, axis=1).values for ct in cts]
        tasks = zip(seeds, samples)
        samples = list(pool.imap_unordered(make_dirichlet, tasks, chunksize=10))
        if close_pool is True:
            pool.close()
        return self.dpca.transform(np.concatenate(samples))

    def fit_transformer(self, iterations=None, samps_per_iteration=None):
        if samps_per_iteration is None:
            samps_per_iteration = 10*self.n_components
        if iterations is None:
            iterations = np.ceil(self.total_samples/samps_per_iteration).astype(int)
        self.pool = mp.Pool()
        for i in tqdm(range(iterations), desc='incremental char pca'):
            samples = self.make_sample(samps_per_iteration,
                                       seed=i+self.seed, pool=self.pool)
            self.transformer.partial_fit(samples)
        del self.transformed_expr
        self.pool.close()
        del self.pool
        self.is_fitted = True
        return self.transformer


class Characterizer2(Characterizer):

    def __init__(self, dirichlet_pca, expressions=None,
                 maxtime=60*60*24*2, maxiter=int(1e5), on_gpu=None,
                 verbose_obtimizer=False, x0=None):
        assert isinstance(dirichlet_pca, Dirichlet_PCA), 'dirichlet_pca musst be of class Dirichlet_PCA'
        self.dpca = dirichlet_pca
        if expressions is None:
            self.expr = dirichlet_pca.expr
            self._new_expr = False
        else:
            assert isinstance(expressions, Expressions), 'experssions musst be of class Expressions'
            self.expr = expressions
            self._new_expr = True
        self.maxtime = maxtime
        self.maxiter = maxiter
        if on_gpu is not None:
            assert isinstance(on_gpu, bool), 'on_gpu musst be boolean'
            self.on_gpu = on_gpu
        self.verbose_obtimizer = verbose_obtimizer
        self.x0 = x0

    @property_deps('cost_dtype', 'devices')
    def on_gpu(self):
        return tf.test.is_gpu_available()

    @property_deps('dimensions')
    def cost_dtype(self):
        if self.on_gpu is True:
            return np.float32
        else:
            return np.float64

    @property_deps()
    def devices(self):
        if self.on_gpu:
            return ['/device:CPU:0', '/device:GPU:0', '/device:GPU:1']
        else:
            return ['/device:CPU:0', '/device:CPU:0', '/device:CPU:0']

    @property_deps('dimensions')
    def x0(self):
        return None

    @property_deps('manifold')
    def dimensions(self):
        if self.x0 is None:
            return 10
        else:
            return self.x0.shape[1]

    @property_deps('transformer')
    def manifold(self):
        return Stiefel(self.dpca.transformer.n_components, self.dimensions)

    @property_deps('transformer')
    def indicator(self):
        data = self.dpca.dataframe['cell type'].values
        blacklisted = self.blacklisted
        indic = dict()
        for ct in self.selection:
            indic[ct] = np.zeros(len(blacklisted))
            ind = (data == ct) & (~blacklisted)
            indic[ct][ind] = 1
            ind = (data != ct) & (~blacklisted)
            indic[ct][ind] = -1
        return indic

    @property_deps()
    def pre_chars(self):
        data = self.dpca.dataframe
        transformed_expr = self.dpca.transformed_expr
        blacklisted = self.blacklisted
        def fit_normal(value, field='cell type'):
            ind = (data[field].values == value) & (~blacklisted)
            dat = transformed_expr[ind, :]
            mean = np.mean(dat, axis=0)
            sigma = np.cov(dat, rowvar=False)
            return value, mean, sigma
        sel = self.selection
        chars = dict()
        with mp.Pool(4) as pool:
            for key, mean, sigma in tqdm(pool.imap_unordered(fit_normal, sel),
                                         total=len(sel), desc='pre character'):
                chars[key] = dict({'mean':mean, 'sigma':sigma})
        return chars

    def cost(self, A):
        data = self.dpca.dataframe
        transformed_expr = self.dpca.transformed_expr
        indic = self.indicator
        t = self.cost_dtype
        #ll = tf.constant(0.0, dtype=t)
        ll = list()
        for ct in self.pre_chars.keys():
            mean = self.pre_chars[ct]['mean'].astype(t)
            sigma = self.pre_chars[ct]['sigma'].astype(t)
            delta = transformed_expr.astype(t) - mean
            with tf.device(self.devices[1]):
                sigprime = tf.linalg.matmul(tf.linalg.matmul(A, sigma, transpose_a=True), A)
                const = tf.linalg.logdet(sigprime)
            with tf.device(self.devices[0]):
                xprime = tf.linalg.matmul(delta, A)
            with tf.device(self.devices[2]):
                var = tf.math.reduce_sum(tf.multiply(tf.linalg.matmul(xprime, sigprime), xprime), axis=1)
            with tf.device(self.devices[0]):
                pre_cost = tf.add(var, const)
                ll.append(tf.math.reduce_logsumexp(tf.math.multiply(pre_cost, indic[ct].astype(t))))
        with tf.device(self.devices[0]):
            result = tf.math.reduce_logsumexp(ll)
        return result

    @property_deps('transforming_A')
    def transformer(self):
        shape = (self.dpca.transformer.n_components, self.dimensions)
        A = tf.Variable(tf.compat.v1.placeholder(self.cost_dtype, shape=shape))
        cost = self.cost(A)
        return pman.Problem(manifold=self.manifold, cost=cost, arg=A)

    @property_deps('transformed_expr', 'affine_trans')
    def transforming_A(self):
        with SuppressPrints(not self.verbose_obtimizer):
            result = SteepestDescent(maxtime=self.maxtime,
                                    maxiter=self.maxiter).solve(self.transformer,
                                                               x=self.x0)
        return result

    def transform(self, values):
        return np.dot(values, self.transforming_A)

    @property_deps()
    def affine_trans(self):
        if not hasattr(self.dpca, 'is_fitted') or self.dpca.is_fitted is not True:
            self.dpca.fit_transformer()
        t1 = self.dpca.transformer
        b1 = t1.mean_
        A1inv = t1.components_
        A1 = A1inv.T
        if t1.whiten:
            A1 = np.dot(A1, np.diag(1/np.sqrt(t1.explained_variance_)))
            A1inv = np.dot(np.diag(np.sqrt(t1.explained_variance_)), A1inv)
        A2 = self.transforming_A
        A2inv = A2.T
        b = b1
        A = np.dot(A1, A2)
        Ainv = np.dot(A2inv, A1inv)
        return b, A, Ainv


class Characterizer3(Characterizer):

    @property_deps('transformed_expr', 'is_fitted', 'affine_trans')
    def transformer(self):
        base = np.zeros((1, self.expr.tau_projected_values[0].shape[1]))
        b = self.dpca.transform(base)
        A = self.pre_means - b
        return b, A.T

    def fit_transformer(self):
        return self.transformer

    def transform(self, values):
        b, A = self.transformer
        return np.dot(values-b, A)

    @property_deps()
    def affine_trans(self):
        if not hasattr(self.dpca, 'is_fitted') or self.dpca.is_fitted is not True:
            self.dpca.fit_transformer()
        t1 = self.dpca.transformer
        b1 = t1.mean_
        A1inv = t1.components_
        A1 = A1inv.T
        if t1.whiten:
            A1 = np.dot(A1, np.diag(1/np.sqrt(t1.explained_variance_)))
            A1inv = np.dot(np.diag(np.sqrt(t1.explained_variance_)), A1inv)
        b2, A2 = self.transformer
        b = b1 + np.dot(b2, A1inv)
        A = np.dot(A1, A2)
        A2inv = np.linalg.pinv(A2)
        Ainv = np.dot(A2inv, A1inv)
        return b, A, Ainv


class Characterizer4(Characterizer):

    @property_deps('transformed_expr', 'is_fitted', 'affine_trans')
    def transformer(self):
        center = np.mean(self.pre_means, axis=0)
        star = self.pre_means.T - center[:, None]
        ray_lengths = np.linalg.norm(star, axis=0)
        center_p = center
        for i in range(self.pre_means.shape[0]-1):
            center_p -= np.dot(center_p, self.pre_means[i,:].T)
        center_p /= np.linalg.norm(center_p) / np.mean(ray_lengths)
        b = center - center_p
        A = self.pre_means.T - b[:, None]
        return b, A

    def fit_transformer(self):
        return self.transformer

    def transform(self, values):
        b, A = self.transformer
        return np.dot(values-b, A)

    @property_deps()
    def affine_trans(self):
        if not hasattr(self.dpca, 'is_fitted') or self.dpca.is_fitted is not True:
            self.dpca.fit_transformer()
        t1 = self.dpca.transformer
        b1 = t1.mean_
        A1inv = t1.components_
        A1 = A1inv.T
        if t1.whiten:
            A1 = np.dot(A1, np.diag(1/np.sqrt(t1.explained_variance_)))
            A1inv = np.dot(np.diag(np.sqrt(t1.explained_variance_)), A1inv)
        b2, A2 = self.transformer
        b = b1 + np.dot(b2, A1inv)
        A = np.dot(A1, A2)
        A2inv = np.linalg.pinv(A2)
        Ainv = np.dot(A2inv, A1inv)
        return b, A, Ainv


def save(data, out):
    with open(out, 'wb') as buff:
        pickle.dump(data, buff, protocol=4)
        
bad_gses = [
    'GSE79706', # most monocytes seem like macrophages
    'GSE41914', # only cohort wrongly classified as naive B-cells while labled as NK-cells
    'GSE116672', # viscRNA-Seq that are very inconsistent
    'GSE66360', # by far most wrongly classified as CD4 T-cells while labled as endothelia cells
    'GSE31312' # by far most wrongly classified as endothelial cells while labled as B-cells
]

if __name__ == "__main__":
    expr = Expressions(grain = args.grain, TSV_dir=args.inFolder)

    if args.center:
        dp = Dirichlet_PCA_centered(expr, samples=args.nSampls,
                                    n_components=args.dcnComp)
    else:
        dp = Dirichlet_PCA(expr, samples=args.nSampls,
                                    n_components=args.dcnComp)
    save({'expr':expr, 'dp':dp}, args.outFile)
    if args.simp:
        ch = Characterizer(dp)
    else:
        ch = Characterizer_centered(dp, n_components=args.nComp)
    bl_ind = ch.data['aggregated_quality'] < args.minquality
    ch.blacklisted = bl_ind
    if args.nonanquality:
        bl_ind |= ch.data['aggregated_quality'].isna()
    #for ct, df in expr.alphas.items():
    #    sample_names = df.columns
    #    nsamp = sum(~ch.blacklisted[sample_names])
    #    print(f'{ct}: {nsamp}')
    bl_ind |= ch.data['gse'].isin(bad_gses)
    ch.blacklisted = bl_ind
    total_bl = sum(ch.blacklisted)
    print(f'{total_bl} samples blacklisted.', file=sys.stderr)

    ch.calculate_probs()
    save({'expr':expr, 'dp':dp, 'ch':ch}, args.outFile)
