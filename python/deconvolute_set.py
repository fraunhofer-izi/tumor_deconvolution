# -*- coding: utf-8 -*-
#
#    Copyright © 2019 Gesellschaft zur Förderung der angewandten Forschung e.V.
#    acting on behalf of its Fraunhofer Institut für Zelltherapie und Immunologie.
#    All rights reserved. Contact: dominik.otto@izi.fraunhofer.de
#
#    This program is free software; you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the Free
#    Software Foundation; either version 3 of the License, or (at your option)
#    any later version.
#
#    This program is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, see <http://www.gnu.org/licenses/>.

import os
import sys
import warnings
import time

# +
tl_vars = ['END', 'CORES', 'sample_count_file']
if all([v in os.environ for v in tl_vars]):
    limited = True
    time_started = time.time()
    running_threads = int(os.environ['CORES'])
    time_end = int(os.environ['END'])
    count_file = os.environ['sample_count_file']
    finish_time_buffer = 600
else:
    limited = False

if 'SLOT' in os.environ:
    slot = os.environ['SLOT']
else:
    slot = 'single'
if 'MODE' in os.environ:
    mode = os.environ['MODE']
else:
    mode = 'normal'
os.environ['OMP_NUM_THREADS'] = '1' # single thread openMP and BLASimport re
os.environ['THEANO_FLAGS'] = 'base_compiledir=~/theano/slot_' + str(slot)
# -

import numpy as np
import pickle
from scipy.stats import multivariate_normal
import pymc3 as pm
import pandas as pd
import theano
import theano.tensor as tt


def load_chars(file):
    with open(file, 'rb') as buff:
        data = pickle.load(buff)
    # should hold: features, A, Ainv, b, chars, gene_means, gene_sd
    return data


def load_sample(csv_file, scale):
    expr = pd.read_csv(csv_file)
    expr = expr.set_index(expr.columns[0])
    sample = sample[np.all(np.isfinite(sample), axis=0),]
    if scale == 'Log':
        sample = np.exp(sample)
    if scale == 'Log2':
        sample = np.exp2(sample)
    if scale == 'Log10':
        sample = 10**sample
    if scale != 'Linear' and scale != 'linear':
        sample -= np.min(sample)
    return sample

def get_grains(file='grains'):
    with open(file) as f:
        grains = list(f)
    return [g.strip() for g in grains]


def decomp_from_trace(trace, old_samps=None):
    if old_samps is None:
        samps = trace['decomp']
    else:
        samps = np.concatenate([old_samps, trace['decomp']], axis=0)
    return np.mean(samps, axis=0), samps


def tau_inv(log_draw):
    o_samp = log_draw - (log_draw.sum() / len(log_draw))
    return o_samp


class CheckAndSave():

    def __init__(self, every=100, tolerance=1e-3,
                 ord=np.inf, save_every=50,
                 save_function=None):
        self.ord = ord
        self.every = every
        self.prev = None
        self.tolerance = tolerance
        self.save_every = save_every
        self.save_i = 0
        if limited is True:
            self.remaining_time = time_end - time_started
        if save_function is None:
            self.save_function = self._pass
        else:
            self.save_function = save_function

    def _pass(*args):
        pass

    def __call__(self, approx, _, i):
        if i % self.every:
            return
        if self.prev is None:
            self.prev = self.get_decomp(approx)
            return
        current = self.get_decomp(approx)
        prev = self.prev
        delta = np.abs(current - prev)
        self.prev = current
        norm = np.linalg.norm(delta, self.ord)
        if norm < self.tolerance:
            self.save_function(self.get_decomp(approx))
            raise StopIteration('Convergence achieved at %d' % i)
        if limited is True:
            try:
                with open(count_file, 'r') as f:
                    remaining_samples = int(f.read())
            except:
                return
            rem_samp_per_thread = np.ceil(remaining_samples / running_threads)
            time_per_samp = self.remaining_time / rem_samp_per_thread
            running = time.time() - time_started + finish_time_buffer
            if running > time_per_samp:
                self.save_function(self.get_decomp(approx))
                raise StopIteration('Time limit reached at %d' % i)
        if not self.save_i % self.save_every:
            self.save_function(self.get_decomp(approx))
        self.save_i += 1

    @staticmethod
    def get_decomp(approx):
        return approx.bij.rmap(approx.mean.get_value())['decomp_stickbreaking__']


class outMeta:

    def __init__(self, set_name, sample_name, out_file, grains=None):
        if grains is None:
            grains = get_grains()
        self.grains = grains
        self.out_file = out_file
        self.sample_names = sample_names
        self.set_name = set_name


class saveFunction:

    def __init__(self, meta, cell_types, grains=None):
        assert isinstance(meta, outMeta), '`meta` must be of class outMeta.'
        self.m = meta
        self.tf = pm.distributions.transforms.StickBreaking()
        self.cell_types = cell_types
        n = len(self.m.grains)
        m = len(self.m.sample_names)
        self.df = pd.DataFrame({'dataset.name':[self.m.set_name]*n*m,
                                'sample.id':[self.m.sample_names]*n,
                                'cell.type':self.m.grains,
                                'prediction':[0.0]*n})
        self.df.index = self.m.grains
        self.index = [ct in self.m.grains for ct in self.cell_types]
        self.out_types = [ct for ct in self.cell_types if ct in self.m.grains]

    def __call__(self, decomp):
        mean_decomp = self.tf.backward(decomp).eval()
        self.df.loc[self.out_types, 'prediction'] = mean_decomp[self.index]
        self.df.to_csv(self.m.out_file, index=False, header=False)


def deconvolute(data, expr, meta=None, sample_deviation=False,
                c_type=None, norm=None, progress=True, max_iter=1e7):
    chars = data['chars']
    features = data['features']
    var = data['variations']
    is_multico = 'multico' in data
    cell_types = list(chars.keys())
    if 'backgrounds' in data:
        background = data['backgrounds']['all']
        for key in data['backgrounds'].keys():
            if key.upper() == c_type.upper():
                background = data[key]
                break
        cell_types += ['background']
    n_types = len(cell_types)
    if meta is not None:
        save_function = saveFunction(meta, cell_types)
    else:
        save_function = None
    test_features = list(expr.index)
    index = np.array([f in test_features for f in features])
    n_features = sum(index)
    filtered_features = [i for (i, v) in zip(features, index) if v]
    sample_float = expr.loc[filtered_features].values
    sample = sample_float.astype(int)
    seq_depth = np.sum(sample)
    ct_prior = np.ones(n_types)
    ct_start = ct_prior/np.sum(ct_prior)
    for i, ct in enumerate(cell_types):
        if ct == 'background':
            ct_prior[i] = 1e-1
            ct_start[i] = 1-ct_start[i]
            ct_start = ct_start/np.sum(ct_start)
            break
    if mode == 'debug':
        mess = '... using {} of {} features with {} available ...'
        print(mess.format(n_features, str(len(index)), str(len(test_features))))
        print('{} reads in total.'.format(sum(expr)))

    def combine_and_embed(samp, deviation, A, Ainv, b):
        base = samp - tt.dot(deviation, A)
        return tt.dot(base, Ainv) + deviation + b
    
    def mix(components, decomp):
        return tt.dot(decomp[None, :], tt.nnet.softmax(components))

    def reduce(samp, A, b):
        return np.dot(samp-b, A)

    def embedd(samp, Ainv, b):
        return np.dot(samp, Ainv) + b

    def project(deviation, A, Ainv, b):
        mapping = reduce(deviation, A, b)
        co_projection = deviation - embedd(mapping, Ainv, b)
        return mapping, co_projection

    l_alphas = sample_float + 1
    t_samp = tau_inv(np.log(l_alphas) - np.log(l_alphas.sum()))
    if is_multico is False:
        A = data['A'][index,:]
        Ainv = data['Ainv'][:,index]
        b = data['b'][index]
        s_Ainv = theano.shared(Ainv)
        s_A = theano.shared(A)
        s_b = theano.shared(b)
        samp_mapping, dev_start = project(t_samp, A, Ainv, b)
    del data

    s = 1
    with pm.Model() as model:
        decomp = pm.Dirichlet('decomp', ct_prior, shape=ct_prior.shape, testval=ct_start)
        cert = (1/(1-decomp) -1)**2
        ct_expr = list()
        scale = pm.Lognormal('scale', testval=10)
        i = 0
        for cell_type in cell_types:
            if cell_type == 'background':
                dev_samp = pm.Normal('comb '+cell_type, mu=background['mean'][index],
                                     sigma=background['std'][index]/scale, shape=(1, n_features),
                                     testval=t_samp)
                ct_expr.append(dev_samp)
                continue
            if is_multico is True:
                A = chars[cell_type]['A'][index,:]
                Ainv = chars[cell_type]['Ainv'][:,index]
                b = chars[cell_type]['b'][index]
                s_Ainv = theano.shared(Ainv)
                s_A = theano.shared(A)
                s_b = theano.shared(b)
                samp_mapping, dev_start = project(dev_start, A, Ainv, b)
                n = A.shape[1]
                samp = pm.Normal(cell_type, sigma=scale/cert[i], shape=(1, n))
            else:
                samp = pm.MvNormal(cell_type, chars[cell_type]['mean'],
                                   cov=chars[cell_type]['sigma']*scale/cert[i],
                                   shape=(1, A.shape[1]),
                                   testval=chars[cell_type]['mean'])
            if sample_deviation is True:
                deviation = pm.Normal('deviation '+cell_type, mu=var[cell_type]['mean'][index],
                                      sigma=var[cell_type]['std'][index]*scale/cert[i],
                                      shape=(1, n_features),
                                      testval=dev_start)
            else:
                deviation = theano.shared(dev_start)
            dev_samp = pm.Deterministic('comb '+cell_type,
                                        combine_and_embed(samp, deviation, s_A, s_Ainv, s_b))
            ct_expr.append(dev_samp)
            i += 1
        ct_expr = tt.concatenate(ct_expr, axis=0)
        transcriptome = pm.Deterministic('trans', mix(ct_expr, decomp))
        obs = pm.Multinomial('obs', seq_depth, transcriptome, observed=sample)
        mean_field = pm.fit(method='advi', n=int(max_iter), progressbar=progress,
                           callbacks=[CheckAndSave(save_function=save_function)],
                           obj_optimizer=pm.adam())

    vals = mean_field.bij.rmap(mean_field.mean.get_value())
    if 'scale_log__' in vals:
        print('Scale {}'.format(np.exp(vals['scale_log__'])))
    decomp = vals['decomp_stickbreaking__']
    return decomp, mean_field


def save_approx(approx, outFile):
    with open(outFile, 'wb') as buff:
        pickle.dump(approx, buff)


def self_test_docker():
    char_file = 'chars.pkl'
    set_name = 'ds1'
    sample_name = 'Sample_2'
    set_file = '/fast_lane_dir/ds1_ensg.csv'
    scale = 'Linear'
    out_file = '/dev/null'
    expr = load_sample(set_file, sample_name, scale)
    data = load_chars(char_file)
    meta = outMeta(set_name, sample_name, out_file)
    decomp, approx = deconvolute(data, expr, meta=meta)
    save_approx(approx, out_file)


def self_test_local():
    char_file = '../docker/data/chars_coarse.pkl'
    set_name = 'ds1'
    sample_name = 'Sample_2'
    set_file = '../docker/data/fast_lane_dir/ds1.csv'
    scale = 'Linear'
    out_file = 'test.csv'
    approx_file = 'approx.pkl'
    grains = get_grains('../docker/data/coarse_grains')
    expr = load_sample(set_file, sample_name, scale)
    data = load_chars(char_file)
    meta = outMeta(set_name, sample_name, out_file, grains=grains)
    decomp, approx = deconvolute(data, expr, meta=meta)
    save_approx(approx, approx_file)


def pre_compile():
    char_file = 'chars.pkl'
    sample_name = 'Sample_1'
    set_file = 'pre.csv'
    scale = 'Linear'
    expr = load_sample(set_file, sample_n2ame, scale)
    data = load_chars(char_file)
    decomp, approx = deconvolute(data, expr, max_iter=1)


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        self_test_docker()
        exit(0)
    if sys.argv[1] == 'precompile':
        print('Precompiling slot {} for mode {}...'.format(slot, mode))
        pre_compile()
        exit(0)
    warnings.filterwarnings('ignore', category=FutureWarning)
    char_file = sys.argv[1]
    set_name = sys.argv[2]
    set_file = sys.argv[3]
    sample_name = sys.argv[4]
    scale = sys.argv[5]
    norm = sys.argv[6]
    c_type = sys.argv[7]
    out_file = sys.argv[8]
    direct_out = sys.argv[9]
    if slot == 'single':
        progress = True
    else:
        progress = False
    if mode == 'debug':
        print('Starting {} of {} ...'.format(sample_name, set_name))
    expr = load_sample(set_file, sample_name, scale)
    data = load_chars(char_file)
    meta = outMeta(set_name, sample_name, out_file)
    decomp, approx = deconvolute(data, expr, meta=meta,
                                 c_type=c_type, norm=norm, progress=progress)
    if mode == 'debug':
        approx_file = '{}/approx_{}_{}.pkl'.format(direct_out, set_name, sample_name)
        save_approx(approx, approx_file)


