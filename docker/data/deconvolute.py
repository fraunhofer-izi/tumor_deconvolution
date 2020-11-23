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
from scipy.special import softmax
import pandas as pd
import theano
theano.config.compute_test_value = 'off'
import theano.tensor as tt
from stickbreaking import *
import pymc3 as pm


def load_chars(file):
    with open(file, 'rb') as buff:
        data = pickle.load(buff)
    return data


def load_sample(csv_file, sample_name, scale, platform=None,
                norm=None, array_scale=5e7, PM_rescale=None):
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
    if (
            norm is not None and
            PM_rescale is not None and
            ((norm.upper() == 'CPM') or (norm.upper() == 'TPM'))
    ):
        total = np.sum(sample)
        if mode == 'debug':
            print('{} reads in total rescaling to {}.'.format(total, PM_rescale))
        if PM_rescale is not None:
            sample *= PM_rescale / total
    elif (platform is not None and
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


def tau_inv_tt(log_draw):
    o_samp = log_draw - (tt.sum(log_draw, axis=1, keepdims=True) / log_draw.shape[1])
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
        rem_samp_per_thread = np.ceil(.5 + (remaining_samples/running_threads))
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
        
    def save(self, approx, sample=False):
        self.save_function(self.get_decomp(approx, sample=sample))

    @staticmethod
    def get_decomp(approx, sample=False):
        if isinstance(approx, pm.variational.approximations.Empirical):
            ap = approx.params[0].eval()
            ap = np.mean(ap, axis=0)
            sb_decomb = approx.bij.rmap(ap)['decomp_stickbreaking__']
        elif sample or isinstance(approx, pm.variational.approximations.NormalizingFlow):
            sb_decomb = approx.sample(10000)['decomp_stickbreaking__']
        elif isinstance(approx, pm.variational.approximations.MeanField):
            ap = approx.mean.get_value()
            sb_decomb = approx.bij.rmap(ap)['decomp_stickbreaking__']
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
        elif len(decomp) == len(self.index):
            mean_decomp = decomp
        else:
            mean_decomp = self.tf.backward(decomp).eval()
        self.df.loc[self.out_types, 'prediction'] = mean_decomp[self.index]
        self.df.to_csv(self.m.out_file, index=False, header=False)


def largest_contributor(A):
    ''' Returns overall unique index of the largest value per column. '''
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


def deconvolute(data, expr, meta=None, sample_deviation=False, resample=False,
                samp_scale=False, method='advi', obj_n_mc=1, nf_obj_n_mc=100,
                svgd_kwargs=dict(n_particles=1000), n_increments=3, project_span=False, poject_samp=True,
                constrain_mix=False, c_type=None, progress=True, init_iter=2e4, max_iter=2e4,
                nfvi_formula='scale-loc-radial*8', sample_background=False):
    assert sample_deviation is False or poject_samp is False
    chars = data['chars']
    features = data['features']
    var = data['variations']
    is_multico = 'multico' in data
    assert is_multico is False or poject_samp is False
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
    #rindex = np.array([f in features for f in test_features])
    n_features = sum(index)
    filtered_features = [i for (i, v) in zip(features, index) if v]
    sample_float = expr.loc[filtered_features].values
    #sample = sample_float.astype(int)
    #seq_depth = np.sum(sample)
    ct_prior = np.ones(n_types)
    ct_start = ct_prior/np.sum(ct_prior)
    for i, ct in enumerate(cell_types):
        if ct == 'background':
            ct_prior[i] = .5
            ct_start[i] = 1-ct_start[i]
            ct_start = ct_start/np.sum(ct_start)
            break
    #save_function.tf = StickBreaking5(ct_prior.shape[0])
    mess = '... using {} of {} features with {} available ...'
    print(mess.format(n_features, str(len(index)), str(len(test_features))))
    n_features = len(index)

    def combine_and_embed(samp, deviation, A, Ainv, b, index):
        base = samp - tt.dot(deviation, A)
        shift = deviation[:, index] + b[:, index]
        result = tt.dot(base, Ainv)[:, index] + shift
        return tau_inv_tt(result)

    def mix(components, decomp):
        return tt.dot(decomp[None, :], tt.nnet.softmax(components))
        
    def reduce(samp, A, b):
        return np.dot(samp-b, A)

    def embedd(samp, Ainv, b):
        return np.dot(samp, Ainv) + b

    def project(sample, A, Ainv, b):
        mapping = reduce(sample, A, b)
        co_projection = sample - embedd(mapping, Ainv, b)
        return mapping, co_projection[None, :]

    l_alphas = sample_float + 1
    magni = l_alphas.sum()
    de_samp = tau_inv(np.log(l_alphas) - np.log(magni))
   
    def make_model(scale=1, dims=slice(None), prior=1, dev_slack=1):
        if is_multico is False:
            A = data['A'][:,dims]
            Ainv = data['Ainv'][dims,:]
            b = data['b']
            s_Ainv = theano.shared(Ainv)
            s_A = theano.shared(A)
            s_b = theano.shared(b[None, :])
            t_samp = np.copy(b)
            t_samp[index] = de_samp
            if project_span:
                a = np.stack([softmax(embedd(ct['mean'], Ainv, b))[index] for ct in chars.values()])
                Q, R = np.linalg.qr(a.T)
                p_samp = Q.T.dot(l_alphas/magni).dot(Q.T)
                dist = pm.Dirichlet.dist(p_samp+1)
                dev_start = np.zeros(b[None, :].shape)
            elif poject_samp is True and not constrain_mix:
                samp_mapping = reduce(t_samp, A, b)
                t_samp = embedd(samp_mapping , Ainv, b)
                dev_start = np.zeros(b[None, :].shape)
                dist = pm.Dirichlet.dist(softmax(t_samp[index])*magni)
            elif poject_samp is True and constrain_mix:
                em_A = np.exp(Ainv.T) / np.sum(np.exp(Ainv.T), axis=0, keepdims=True)
                em_A_inv = np.linalg.pinv(em_A)
                t_samp = np.dot(np.dot(t_samp, em_A), em_A_inv)
                dev_start = np.zeros(b[None, :].shape)
                dist = pm.Dirichlet.dist(t_samp[index]*magni)
            else:
                samp_mapping, dev_start = project(t_samp, A, Ainv, b)
                dev_start = np.zeros(b[None, :].shape)
                dist = pm.Dirichlet.dist(l_alphas)
                
        with pm.Model() as model:
            decomp = pm.Dirichlet('decomp', ct_prior*prior, shape=ct_prior.shape,
                                  testval=ct_start, transform=save_function.tf)
            ct_expr = list()
            if samp_scale is True:
                scale = pm.Lognormal('scale', testval=10)
            for i, cell_type in enumerate(cell_types):
                if cell_type == 'background':
                    dev_samp = pm.Normal('comb '+cell_type, mu=background['mean'],
                                         sigma=background['std']*dev_slack,
                                         shape=(1, n_features), testval=t_samp)
                    ct_expr.append(dev_samp)
                    continue
                if is_multico is True:
                    A = chars[cell_type]['A'][:,dims]
                    Ainv = chars[cell_type]['Ainv'][dims,:]
                    b = chars[cell_type]['b']
                    s_Ainv = theano.shared(Ainv)
                    s_A = theano.shared(A)
                    s_b = theano.shared(b[None, :])
                    samp_mapping, dev_start = project(t_samp, A, Ainv, b)
                    n = A.shape[1]
                    samp = pm.Normal(cell_type, sigma=scale, shape=(1, n))
                else:
                    samp = pm.MvNormal(cell_type, chars[cell_type]['mean'][dims],
                                       cov=chars[cell_type]['sigma'][dims, dims]*scale,
                                       shape=(1, A.shape[1]),
                                       testval=chars[cell_type]['mean'][dims])
                if sample_deviation is True:
                    deviation = pm.Normal('deviation '+cell_type, mu=dev_start,
                                          sigma=var[cell_type]['std']*dev_slack,
                                          shape=(1, n_features))
                else:
                    deviation = theano.shared(dev_start)
                dev_samp = pm.Deterministic('comb '+cell_type,
                                            combine_and_embed(samp, deviation, s_A, s_Ainv, s_b, index))
                ct_expr.append(dev_samp)
            ct_expr = tt.concatenate(ct_expr, axis=0)
            if constrain_mix and dims == slice(None):
                # project the mix to the subspace spanned by softmax(Ainv.T)
                #em_A = (np.exp(Ainv) / np.sum(np.exp(Ainv), axis=0, keepdims=True)).T
                em_A = tt.nnet.softmax(s_Ainv).T[index,:]
                em_A_inv = tt.nlinalg.MatrixPinv()(em_A)
                def mix(components, decomp):
                    mixed = tt.dot(decomp[None, :], tt.nnet.softmax(components))
                    return tt.dot(tt.dot(mixed, em_A), em_A_inv)
            else:
                def mix(components, decomp):
                    return tt.dot(decomp[None, :], tt.nnet.softmax(components))
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
        if resample:
            decomps = approx.sample(resample)['decomp_stickbreaking__']
            decomp = np.mean(decomps, axis=0)
        else:
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
        if resample:
            decomps = approx.sample(resample)['decomp_stickbreaking__']
            decomp = np.mean(decomps, axis=0)
        else:
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
            if not sf.musst_stop(buffer=300):
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
        if resample:
            decomps = approx.sample(resample)['decomp_stickbreaking__']
            decomp = np.mean(decomps, axis=0)
        else:
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
    elif method == 'map':
        approx = pm.find_MAP(model=model, progressbar=progress)
        decomp = approx['decomp']
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
    set_file = '/input_test/ds1.csv'
    scale = 'Linear'
    out_file = '/dev/null'
    platform = 'AFFYMETRIX'
    expr = load_sample(set_file, sample_name, scale, platform)
    data = load_chars(char_file)
    meta = outMeta(set_name, sample_name, out_file)
    decomp, approx = deconvolute(data, expr, meta=meta,
                                 method='advi', max_iter=10)
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
    expr = load_sample(set_file, sample_name, scale, platform, norm)
    data = load_chars(char_file)
    meta = outMeta(set_name, sample_name, out_file)
    decomp, approx = deconvolute(data, expr, meta=meta,
                                 c_type=c_type, progress=progress)
    if mode == 'debug':
        approx_file = '{}/approx_{}_{}.pkl'.format(direct_out, set_name, sample_name)
        save_approx(approx, approx_file)


