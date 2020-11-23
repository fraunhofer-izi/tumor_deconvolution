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
import time

# +
time_started = time.time()
tl_vars = ['END', 'CORES', 'sample_count_file']
if all([v in os.environ for v in tl_vars]):
    running_threads = int(os.environ['CORES'])
    time_end = int(os.environ['END'])
    count_file = os.environ['sample_count_file']
    finish_time_buffer = 600
    limited = time_started < time_end
else:
    limited = False

if 'SLOT' in os.environ and os.environ['SLOT'].upper() != 'SINGLE':
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

import sys
import re
import warnings
import numpy as np
import pickle
from scipy.stats import multivariate_normal
import pymc3 as pm
import pandas as pd
import theano
import theano.tensor as tt
from stickbreaking import *


def load_chars(file):
    with open(file, 'rb') as buff:
        data = pickle.load(buff)
    return data


def load_sample(csv_file, sample_name, scale, platform=None, array_scale=5e7):
    expr = pd.read_csv(csv_file)
    expr = expr.set_index(expr.columns[0])
    if sample_name not in expr.columns:
        raise ValueError('Samplename not in data columns.')
    sample = expr[sample_name]
    sample = sample[np.isfinite(sample)]
    if scale == 'Log':
        sample = np.exp(sample)
    if scale == 'Log2':
        sample = np.exp2(sample)
    if scale == 'Log10':
        sample = 10**sample
    if (scale != 'Linear' and scale != 'linear') or any(sample<0):
        sample -= np.min(sample)
    if (platform is not None and
            (
                ('AFFYMETRIX' in platform.upper()) or
                ('AGILENT' in platform.upper()) or
                ('OPERON' in platform.upper())
            )
       ):
        total = np.sum(sample)
        if mode == 'debug':
            print('{} reads in total rescaling to {}.'.format(total, array_scale))
        if array_scale is not None:
            sample *= array_scale / total
    elif mode == 'debug':
        print('{} reads in total.'.format(np.sum(sample)))
    return sample

def get_grains(file='grains'):
    with open(file) as f:
        grains = list(f)
    return [g.strip() for g in grains]


def decomp_from_trace(trace, old_samps=None):
    if old_samps is None:
        samps = trace['decomp_stickbreaking__']
    else:
        samps = np.concatenate([old_samps, trace['decomp_stickbreaking__']], axis=0)
    return np.mean(samps, axis=0)


def tau_inv(log_draw):
    o_samp = log_draw - (log_draw.sum() / len(log_draw))
    return o_samp


class CheckAndSave():

    def __init__(self, every=100, tolerance=1e-5,
                 ord=np.inf, save_every=50, save_function=None):
        self.ord = ord
        self.every = every
        self.prev = None
        self.tolerance = tolerance
        self.save_every = save_every
        self.save_i = 0
        self.limited = limited
        self.last_time = None
        if self.limited is True:
            self.remaining_time = time_end - time_started
        if save_function is None:
            self.save_function = self._pass
        else:
            self.save_function = save_function

    def _pass(*args):
        pass
    
    def musst_stop(self, checkpoint=True, buffer=None):
        if self.limited is not True:
            return False
        try:
            with open(count_file, 'r') as f:
                remaining_samples = int(f.read())
        except:
            return False
        rem_samp_per_thread = np.ceil(remaining_samples / running_threads)
        time_per_samp = self.remaining_time / rem_samp_per_thread
        now = time.time()
        running = now - time_started + finish_time_buffer
        if buffer is None and self.last_time is not None:
            since_last = now - self.last_time
            buffer = since_last/2
        elif buffer is None:
            buffer = 0
        if checkpoint is True:
            self.last_time = now
        return running + buffer > time_per_samp

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
            self.save(approx)
            raise StopIteration('Convergence achieved at %d' % i)
        if self.musst_stop():
            self.save(approx)
            raise StopIteration('Time limit reached at %d' % i)
        if not self.save_i % self.save_every:
            self.save(approx)
        self.save_i += 1
        
    def save(self, approx):
        self.save_function(self.get_decomp(approx))

    @staticmethod
    def get_decomp(approx):
        if isinstance(approx, pm.variational.approximations.Empirical):
            ap = approx.params[0].eval()
            ap = np.mean(ap, axis=0)
            sb_decomb = approx.bij.rmap(ap)['decomp_stickbreaking__']
        elif isinstance(approx, pm.variational.approximations.MeanField):
            ap = approx.mean.get_value()
            sb_decomb = approx.bij.rmap(ap)['decomp_stickbreaking__']
        elif isinstance(approx, pm.variational.approximations.NormalizingFlow):
            sb_decomb = approx.sample(10000)['decomp_stickbreaking__']
        else:
            message = 'Save method not implemented for approximations of type {}'
            raise NotImplementedError(message.format(type(approx)))
        if len(sb_decomb.shape) > 1:
            return np.mean(sb_decomb, axis=0)
        return sb_decomb


class outMeta:

    def __init__(self, set_name, sample_name, out_file, grains=None):
        if grains is None:
            grains = get_grains()
        self.grains = grains
        self.out_file = out_file
        self.sample_name = sample_name
        self.set_name = set_name


class saveFunction:

    def __init__(self, meta, cell_types, grains=None):
        assert isinstance(meta, outMeta), '`meta` must be of class outMeta.'
        self.m = meta
        self.tf = pm.distributions.transforms.StickBreaking()
        self.cell_types = cell_types
        n = len(self.m.grains)
        self.df = pd.DataFrame({'dataset.name':[self.m.set_name]*n,
                                'sample.id':[self.m.sample_name]*n,
                                'cell.type':self.m.grains,
                                'prediction':[0.0]*n})
        self.df.index = self.m.grains
        self.index = [ct in self.m.grains for ct in self.cell_types]
        self.out_types = [ct for ct in self.cell_types if ct in self.m.grains]

    def __call__(self, decomp):
        if isinstance(decomp, pm.backends.base.MultiTrace):
            mean_decomp = decomp_from_trace(decomp)
        else:
            mean_decomp = self.tf.backward(decomp).eval()
        self.df.loc[self.out_types, 'prediction'] = mean_decomp[self.index]
        self.df.to_csv(self.m.out_file, index=False, header=False)


def largest_contributor(A):
    ''' Returns unique index of the largest value per column. '''
    largest = np.argmax(A, axis=0)
    while len(set(largest)) < len(largest):
        seen = set()
        dupes = set()
        for val in largest:
            if val in seen:
                dupes.add(val)
            seen.add(val)
        for dupe in dupes:
            replace = largest == dupe
            best_index = np.argmax(A[dupe, replace])
            best_index = np.where(replace)[0][best_index]
            replace[best_index] = False
            mask = np.ones(A.shape[0], np.bool)
            mask[largest[~replace]] = False
            largest[replace] = np.argmax(A[:, replace][mask, :], axis=0)
    return largest


char_file = '../docker/data/chars.pkl'
set_name = 'DS487'
sample_name = 'S93'
set_file = '../docker/input_synapse2/DS487-hugo-gene-expr.csv'
scale = 'Log2'
out_file = 'test.csv'
approx_file = 'approx.pkl'
grains = get_grains('../docker/data/coarse_grains')
expr = load_sample(set_file, sample_name, scale)
data = load_chars(char_file)
meta = outMeta(set_name, sample_name, out_file, grains=grains)

# +
meta=None
sample_deviation=True
samp_scale=False
method='advi'
obj_n_mc=8
nf_obj_n_mc=100
svgd_kwargs=dict(n_particles=4000)
n_increments=3
c_type=None
norm=None
progress=True
init_iter=32e3
max_iter=1e5
nfvi_formula='scale-loc-radial*8'
sample_background=False

chars = data['chars']
features = data['features']
var = data['variations']
is_multico = 'multico' in data
cell_types = list(chars.keys())
if 'backgrounds' in data and sample_background is True:
    background = data['backgrounds']['all']
    if c_type is not None:
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
#sample = sample_float.astype(int)
#seq_depth = np.sum(sample)
ct_prior = np.ones(n_types)
ct_start = ct_prior/np.sum(ct_prior)
for i, ct in enumerate(cell_types):
    if ct == 'background':
        ct_prior[i] = 5e-1
        ct_start[i] = 1-ct_start[i]
        ct_start = ct_start/np.sum(ct_start)
        break
mess = '... using {} of {} features with {} available ...'
print(mess.format(n_features, str(len(index)), str(len(test_features))))

def combine_and_embed(samp, deviation, A, Ainv, b):
    base = samp - tt.dot(deviation, A)
    return tt.dot(base, Ainv) + deviation + b

def mix(components, decomp):
    return tt.dot(decomp[None, :], tt.nnet.softmax(components))

def reduce(samp, A, b):
    return np.dot(samp-b, A)

def embedd(samp, Ainv, b):
    return np.dot(samp, Ainv) + b

def project(sample, A, Ainv, b):
    mapping = reduce(sample, A, b)
    co_projection = sample - embedd(mapping, Ainv, b)
    return mapping, co_projection

l_alphas = sample_float + 1
t_samp = tau_inv(np.log(l_alphas) - np.log(l_alphas.sum()))
dist = pm.Dirichlet.dist(l_alphas)

def make_model(scale=1e-1, dims=slice(None), prior=10, dev_slack=1e7):
    if is_multico is False:
        A = data['A'][index,dims]
        Ainv = data['Ainv'][dims,index]
        b = data['b'][index]
        s_Ainv = theano.shared(Ainv)
        s_A = theano.shared(A)
        s_b = theano.shared(b)
        samp_mapping, dev_start = project(t_samp, A, Ainv, b)

    with pm.Model() as model:
        decomp = pm.Dirichlet('decomp', ct_prior*prior, shape=ct_prior.shape,
                              testval=ct_start)#, transform=StickBreaking5(ct_prior.shape[0]))
        ct_expr = list()
        if samp_scale is True:
            scale = pm.Lognormal('scale', testval=10)
        for i, cell_type in enumerate(cell_types):
            if cell_type == 'background':
                dev_samp = pm.Normal('comb '+cell_type, mu=background['mean'][index],
                                     sigma=background['std'][index]/scale, shape=(1, n_features),
                                     testval=t_samp)
                ct_expr.append(dev_samp)
                continue
            if is_multico is True:
                A = chars[cell_type]['A'][index,dims]
                Ainv = chars[cell_type]['Ainv'][dims,index]
                b = chars[cell_type]['b'][index]
                s_Ainv = theano.shared(Ainv)
                s_A = theano.shared(A)
                s_b = theano.shared(b)
                samp_mapping, dev_start = project(t_samp, A, Ainv, b)
                n = A.shape[1]
                samp = pm.Normal(cell_type, sigma=scale, shape=(1, n))
            else:
                samp = pm.MvNormal(cell_type, chars[cell_type]['mean'][dims],
                                   cov=chars[cell_type]['sigma'][dims, dims]*scale,
                                   shape=(1, A.shape[1]),
                                   testval=chars[cell_type]['mean'][dims])
            if sample_deviation is True:
                deviation = pm.Normal('deviation '+cell_type, mu=var[cell_type]['mean'][index],
                                      sigma=var[cell_type]['std'][index]*dev_slack,
                                      shape=(1, n_features), testval=dev_start)
            else:
                deviation = theano.shared(dev_start)
            dev_samp = pm.Deterministic('comb '+cell_type,
                                        combine_and_embed(samp, deviation, s_A, s_Ainv, s_b))
            ct_expr.append(dev_samp)
        ct_expr = tt.concatenate(ct_expr, axis=0)
        transcriptome = pm.Deterministic('trans', mix(ct_expr, decomp))
        pot = pm.Potential('obs', dist.logp(transcriptome))
        #obs = pm.Multinomial('obs', seq_depth, transcriptome, observed=sample, dtype='int64')
    return model

sf = CheckAndSave(save_function=save_function)
if method != 'increment':
    if mode == 'debug':
        print('Compiling model ...')
    model = make_model()
# -

map_estimate3 = pm.find_MAP(model=model)
map_estimate3['decomp']

pred = pd.DataFrame({'precited':map_estimate3['decomp'], 'cell.types':cell_types}).set_index('cell.types')
gen_gold_path = ('/<censored_path>/dominik.otto/'
                 + 'tumor-deconvolution-dream-challenge/synapse/'
                 + 'gold_standards/lb_coarse_r2.csv')
gg = pd.read_csv(gen_gold_path).set_index(['dataset.name', 'sample.id', 'cell.type'])
mes = gg.loc[(set_name, sample_name)]
mes.join(pred)

map_estimate3 = pm.find_MAP(model=model)
map_estimate3['decomp']

pred = pd.DataFrame({'precited':map_estimate3['decomp'], 'cell.types':cell_types}).set_index('cell.types')
gen_gold_path = ('/<censored_path>/dominik.otto/'
                 + 'tumor-deconvolution-dream-challenge/synapse/'
                 + 'gold_standards/lb_coarse_r2.csv')
gg = pd.read_csv(gen_gold_path).set_index(['dataset.name', 'sample.id', 'cell.type'])
mes = gg.loc[(set_name, sample_name)]
mes.join(pred)

cell_types

approx2_advi = pm.fit(model=model, n=int(1e4))

advi_samp = approx2_advi.sample(draws=int(1e4))

pred = pd.DataFrame({'precited':np.mean(advi_samp['decomp'], axis=0), 'cell.types':cell_types}).set_index('cell.types')
gen_gold_path = ('/<censored_path>/dominik.otto/'
                 + 'tumor-deconvolution-dream-challenge/synapse/'
                 + 'gold_standards/lb_coarse_r2.csv')
gg = pd.read_csv(gen_gold_path).set_index(['dataset.name', 'sample.id', 'cell.type'])
mes = gg.loc[(set_name, sample_name)]
mes.join(pred)

pm.pairplot(advi_samp, var_names='decomp');

approx_ladvi2 = pm.fit(model=model, n=int(1e5))

approx_ladvi2 = pm.fit(model=model, n=int(1e5))

approx_ladvi3 = pm.fit(model=model, n=int(1e5))

ladvi_samp3 = approx_ladvi2.sample(draws=int(1e4))

pred = pd.DataFrame({'precited':np.mean(ladvi_samp3['decomp'], axis=0), 'cell.types':cell_types}).set_index('cell.types')
gen_gold_path = ('/<censored_path>/dominik.otto/'
                 + 'tumor-deconvolution-dream-challenge/synapse/'
                 + 'gold_standards/lb_coarse_r2.csv')
gg = pd.read_csv(gen_gold_path).set_index(['dataset.name', 'sample.id', 'cell.type'])
mes = gg.loc[(set_name, sample_name)]
mes.join(pred)

pm.pairplot(ladvi_samp2, var_names='decomp');

pred = pd.DataFrame({'precited':np.mean(ladvi_samp2['decomp'], axis=0), 'cell.types':cell_types}).set_index('cell.types')
gen_gold_path = ('/<censored_path>/dominik.otto/'
                 + 'tumor-deconvolution-dream-challenge/synapse/'
                 + 'gold_standards/lb_coarse_r2.csv')
gg = pd.read_csv(gen_gold_path).set_index(['dataset.name', 'sample.id', 'cell.type'])
mes = gg.loc[(set_name, sample_name)]
mes.join(pred)

approx3 = pm.sample(model=model, draws=int(1e5), cores=4,
                  init='advi', n_init=int(1e4), chains=16)
with open('../../data/experimental_decomp_4.pkl', 'wb') as buff:
    pickle.dump(approx3, buff, protocol=4)

# + {"active": ""}
# with open('../../data/experimental_decomp_3.pkl', 'rb') as buff:
#     approx3 = pickle.load(buff)
# -

pm.pairplot(approx3, var_names='decomp');

pm.traceplot(approx3, var_names='decomp');

pm.traceplot(approx3, var_names='monocytic.lineage');

# +
import plotnine as pn

def ct_trace(cell_type):
    trace = approx3.get_values(varname=cell_type)
    df = pd.DataFrame(trace[:, 0, :], columns=[f'PC {i+1}' for i in range(trace.shape[2])])
    df['cell type'] = 'trace'
    return df

def ct_means(cell_type):
    posis = pd.DataFrame()
    for c in approx3.chains:
        position = np.mean(approx3.get_values(varname=cell_type, chains=c), axis=0)
        df = pd.DataFrame(position, columns=[f'PC {i+1}' for i in range(position.shape[1])])
        df['chain'] = str(c)
        posis = posis.append(df)
    return posis

def ct_plot(cell_type, ggplot=None, x='PC 1', y='PC 2', mean=True):
    if ggplot is None:
        ggplot = pn.ggplot(pn.aes(x=x, y=y))
    if mean is True:
        posis = ct_means(cell_type)
        pl = (ggplot 
            + pn.geom_point(pn.aes(color='chain'), data=posis, size=8, alpha=.3)
            + pn.geom_point(pn.aes(color='chain'), data=posis))
    else:
        posis = ct_trace(cell_type)
        pl = (ggplot + pn.geom_point(data=posis))
    return pl


# -

ct_plot('monocytic.lineage', x='PC 4', y='PC 2')

ct_plot('monocytic.lineage', x='PC 4', y='PC 2', mean=False)

file_name = '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/transformations_coarse_HUGO_array.pkl'
with open(file_name, 'rb') as buff:
    data = pickle.load(buff)
ch = data['ch']
del data

ct_plot('endothelial.cells', ch.plot(1, 2, ['monocytic.lineage', 'endothelial.cells'])).draw();

ct_plot('monocytic.lineage', ch.plot(4, 2, ['CD4.T.cells', 'monocytic.lineage']), mean=False).draw();

ct_plot('monocytic.lineage', ch.plot(4, 2, ['CD4.T.cells', 'monocytic.lineage'])).draw();

ct_plot('CD4.T.cells', ch.plot(4, 2, ['CD4.T.cells', 'monocytic.lineage'])).draw();

ct_plot('B.cells', ch.plot(4, 2, ['CD4.T.cells', 'monocytic.lineage', 'B.cells'])).draw();

ct_plot('B.cells', ch.plot(4, 2, ['CD4.T.cells', 'monocytic.lineage', 'B.cells']), mean=False).draw();

ct_plot('neutrophils', ch.plot(4, 2, ['CD4.T.cells', 'monocytic.lineage', 'neutrophils']), mean=False).draw();


def deconvolute(data, expr, meta=None, sample_deviation=False,
                samp_scale=False, method='combination', obj_n_mc=8, nf_obj_n_mc=100,
                svgd_kwargs=dict(n_particles=4000), n_increments=3,
                c_type=None, norm=None, progress=True, init_iter=32e3, max_iter=1e5,
                nfvi_formula='scale-loc-radial*8', sample_background=False):
    chars = data['chars']
    features = data['features']
    var = data['variations']
    is_multico = 'multico' in data
    cell_types = list(chars.keys())
    if 'backgrounds' in data and sample_background is True:
        background = data['backgrounds']['all']
        if c_type is not None:
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
    #sample = sample_float.astype(int)
    #seq_depth = np.sum(sample)
    ct_prior = np.ones(n_types)
    ct_start = ct_prior/np.sum(ct_prior)
    for i, ct in enumerate(cell_types):
        if ct == 'background':
            ct_prior[i] = 5e-1
            ct_start[i] = 1-ct_start[i]
            ct_start = ct_start/np.sum(ct_start)
            break
    mess = '... using {} of {} features with {} available ...'
    print(mess.format(n_features, str(len(index)), str(len(test_features))))

    def combine_and_embed(samp, deviation, A, Ainv, b):
        base = samp - tt.dot(deviation, A)
        return tt.dot(base, Ainv) + deviation + b

    def mix(components, decomp):
        return tt.dot(decomp[None, :], tt.nnet.softmax(components))

    def reduce(samp, A, b):
        return np.dot(samp-b, A)

    def embedd(samp, Ainv, b):
        return np.dot(samp, Ainv) + b

    def project(sample, A, Ainv, b):
        mapping = reduce(sample, A, b)
        co_projection = sample - embedd(mapping, Ainv, b)
        return mapping, co_projection

    l_alphas = sample_float + 1
    t_samp = tau_inv(np.log(l_alphas) - np.log(l_alphas.sum()))
    dist = pm.Dirichlet.dist(l_alphas)

    def make_model(scale=1e-3, dims=slice(None), prior=10):
        if is_multico is False:
            A = data['A'][index,dims]
            Ainv = data['Ainv'][dims,index]
            b = data['b'][index]
            s_Ainv = theano.shared(Ainv)
            s_A = theano.shared(A)
            s_b = theano.shared(b)
            samp_mapping, dev_start = project(t_samp, A, Ainv, b)

        with pm.Model() as model:
            decomp = pm.Dirichlet('decomp', ct_prior*prior, shape=ct_prior.shape,
                                  testval=ct_start)#, transform=StickBreaking5(ct_prior.shape[0]))
            ct_expr = list()
            if samp_scale is True:
                scale = pm.Lognormal('scale', testval=10)
            for i, cell_type in enumerate(cell_types):
                if cell_type == 'background':
                    dev_samp = pm.Normal('comb '+cell_type, mu=background['mean'][index],
                                         sigma=background['std'][index]/scale, shape=(1, n_features),
                                         testval=t_samp)
                    ct_expr.append(dev_samp)
                    continue
                if is_multico is True:
                    A = chars[cell_type]['A'][index,dims]
                    Ainv = chars[cell_type]['Ainv'][dims,index]
                    b = chars[cell_type]['b'][index]
                    s_Ainv = theano.shared(Ainv)
                    s_A = theano.shared(A)
                    s_b = theano.shared(b)
                    samp_mapping, dev_start = project(t_samp, A, Ainv, b)
                    n = A.shape[1]
                    samp = pm.Normal(cell_type, sigma=scale, shape=(1, n))
                else:
                    samp = pm.MvNormal(cell_type, chars[cell_type]['mean'][dims],
                                       cov=chars[cell_type]['sigma'][dims, dims]*scale,
                                       shape=(1, A.shape[1]),
                                       testval=chars[cell_type]['mean'][dims])
                if sample_deviation is True:
                    deviation = pm.Normal('deviation '+cell_type, mu=var[cell_type]['mean'][index],
                                          sigma=var[cell_type]['std'][index]*scale,
                                          shape=(1, n_features), testval=dev_start)
                else:
                    deviation = theano.shared(dev_start)
                dev_samp = pm.Deterministic('comb '+cell_type,
                                            combine_and_embed(samp, deviation, s_A, s_Ainv, s_b))
                ct_expr.append(dev_samp)
            ct_expr = tt.concatenate(ct_expr, axis=0)
            transcriptome = pm.Deterministic('trans', mix(ct_expr, decomp))
            pot = pm.Potential('obs', dist.logp(transcriptome))
            #obs = pm.Multinomial('obs', seq_depth, transcriptome, observed=sample, dtype='int64')
        return model

    sf = CheckAndSave(save_function=save_function)
    if method != 'increment':
        if mode == 'debug':
            print('Compiling model ...')
        model = make_model()
        if mode == 'debug':
            print('Starting inference ...')
    if method == 'increment':
        if is_multico is True:
            message = 'The method increment is not implemented for multico characterisations.'
            raise NotImplementedError(message)
        maxdim = np.min([data['A'].shape[1], 50])
        start = np.min([n_types, maxdim])
        steps = np.unique(np.geomspace(start, maxdim, num=n_increments, dtype=int))
        if mode == 'debug':
            print('Doing increments {} ...'.format(steps))
        lastparam = None
        for dims in steps:
            print('Increment {}'.format(dims))
            if mode == 'debug':
                print('Compiling model ...')
            start_compile = time.time()
            model = make_model(dims=slice(dims))
            compile_time = time.time() - start_compile
            if mode == 'debug':
                print('Starting inference ...')
            with model:
                advi = pm.ADVI()
                if lastparam is not None:
                    rmap = advi.approx.groups[0].bij.rmap
                    newpars = {param.name: rmap(param.eval())
                        for param in advi.approx.params}
                    for ct in cell_types:
                        if ct == 'background':
                            continue
                        mus = lastparam['mu'][ct]
                        ind = np.indices(lastparam['mu'][ct].shape, sparse=True)
                        lastparam['mu'][ct] = newpars['mu'][ct]
                        lastparam['mu'][ct][ind] = mus

                        rohs = lastparam['rho'][ct]
                        ind = np.indices(lastparam['rho'][ct].shape, sparse=True)
                        lastparam['rho'][ct] = newpars['rho'][ct]
                        lastparam['rho'][ct][ind] = rohs

                    fmap = advi.approx.groups[0].bij.map
                    advi.approx.params[0].set_value(fmap(lastparam['mu']))
                    advi.approx.params[1].set_value(fmap(lastparam['rho']))
                approx = advi.fit(n=int(max_iter), progressbar=progress,
                                  callbacks=[sf], obj_n_mc=obj_n_mc)
            if sf.musst_stop(buffer=compile_time+60):
                print('Stopping: Not enought time for next increment ...')
                break
            rmap = approx.groups[0].bij.rmap
            lastparam = {param.name: rmap(param.eval())
                for param in approx.params}
            lastdims = dims
        decomp = lastparam['mu']['decomp_stickbreaking__']
    elif method == 'advi':
        approx = pm.fit(model=model, method='advi', n=int(max_iter),
                        progressbar=progress, callbacks=[sf], obj_n_mc=obj_n_mc)
        vals = approx.bij.rmap(approx.mean.get_value())
        if 'scale_log__' in vals:
            print('Scale {}'.format(np.exp(vals['scale_log__'])))
        decomp = sf.get_decomp(approx)
    elif method == 'decrate':
        approx = pm.fit(model=model, method='advi', n=int(max_iter), obj_optimizer=pm.adam(),
                        progressbar=progress, callbacks=[sf], obj_n_mc=obj_n_mc)
        vals = approx.bij.rmap(approx.mean.get_value())
        if 'scale_log__' in vals:
            print('Scale {}'.format(np.exp(vals['scale_log__'])))
        decomp = sf.get_decomp(approx)
    elif method == 'svgd':
        sf.every = 20
        approx = pm.fit(model=model, method='svgd', inf_kwargs=svgd_kwargs,
                        n=int(max_iter), progressbar=progress, callbacks=[sf],
                        obj_n_mc=obj_n_mc)
        vals = approx.params[0].eval()
        vals = np.mean(vals, aixs=0)
        vals = approx.bij.rmap(vals)
        if 'scale_log__' in vals:
            print('Scale {}'.format(np.exp(vals['scale_log__'])))
        decomp = vals['decomp_stickbreaking__']
    elif method == 'nfvi':
        with model:
            nfvi = pm.NFVI(nfvi_formula)
            approx = nfvi.fit(n=int(max_iter), progressbar=progress,
                              callbacks=[sf], obj_n_mc=nf_obj_n_mc)
        decomps = approx.sample(1000)['decomp_stickbreaking__']
        decomp = np.mean(decomps, axis=0)
    elif method == 'nuts':
        approx = pm.sample(model=model, draws=int(max_iter), progressbar=progress,
                          init='advi', n_init=int(init_iter), chains=1)
        decomp = decomp_from_trace(approx)
    elif method == 'combination':
        with model:
            approx = pm.fit(method='advi', n=int(init_iter), progressbar=progress,
                            callbacks=[sf], obj_n_mc=obj_n_mc)
            print('Starting SVGD ...')
            n = svgd_kwargs.get('n_particles', 100)
            rmap = approx.bij.rmap
            svgd = pm.SVGD(**svgd_kwargs)
            fmap = svgd.approx.bij.map
            means = fmap(rmap(approx.mean.get_value()))
            stds = fmap(rmap(approx.std.eval()))
            start = np.random.normal(means, stds, size=(n, len(means)))
            svgd.approx.params[0].set_value(start)
            sf.every = 20
            sf.last_time = None
            approx = svgd.fit(n=int(max_iter), progressbar=progress,
                              callbacks=[sf], obj_n_mc=obj_n_mc)
        decomp = sf.get_decomp(approx)
    elif method == 'nfcomb':
        assert re.match('scale-loc', nfvi_formula), \
            ('The nfvi formula needs to start with `scale-loc` in order '
             + f'to be initiated with advi. Instead it is `{nfvi_formula}`.')
        with model:
            approx = pm.fit(method='advi', n=int(init_iter), progressbar=progress,
                            callbacks=[sf], obj_n_mc=obj_n_mc)
            print('Starting NFVI ...')
            rmap = approx.bij.rmap
            startpars = {param.name: rmap(param.eval())
                         for param in approx.params}
            nfvi = pm.NFVI(nfvi_formula)
            fmap = nfvi.approx.bij.map
            nfvi.approx.params[-1].set_value(fmap(startpars['rho']))
            nfvi.approx.params[-2].set_value(fmap(startpars['mu']))
            approx = nfvi.fit(n=int(max_iter), progressbar=progress,
                              callbacks=[sf], obj_n_mc=nf_obj_n_mc)
        decomp = sf.get_decomp(approx)
    else:
        message = 'The method {} is not implemented.'
        raise NotImplementedError(message.format(method))
    if save_function is not None:
        save_function(decomp)
    return decomp, approx


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
    char_file = '../docker/data/chars.pkl'
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
    platform = sys.argv[8]
    out_file = sys.argv[9]
    direct_out = sys.argv[10]
    print('Writing to {} ...'.format(direct_out))
    if slot == 'single':
        progress = True
    else:
        progress = False
    print('Starting {} of {} ...'.format(sample_name, set_name))
    expr = load_sample(set_file, sample_name, scale, platform)
    data = load_chars(char_file)
    meta = outMeta(set_name, sample_name, out_file)
    decomp, approx = deconvolute(data, expr, meta=meta,
                                 c_type=c_type, norm=norm, progress=progress)
    if mode == 'debug':
        approx_file = '{}/approx_{}_{}.pkl'.format(direct_out, set_name, sample_name)
        save_approx(approx, approx_file)


