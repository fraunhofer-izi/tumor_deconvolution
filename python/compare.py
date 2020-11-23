import pandas as pd
import _pickle as pickle
import plotnine as pn
import re
import numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.ndimage.filters import gaussian_filter1d
import warnings
import matplotlib.cbook
import os
from tqdm.auto import tqdm
import multiprocess as mp

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=matplotlib.cbook.MatplotlibDeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

class Compare:
    
    gs_paths = {'coarse':[
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/examples/example_files/example_gold_standard/fast_lane_course.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/generated_goldstandard/gs_coarse2.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/generated_goldstandard/gs_hugo3.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/synapse/gold_standards/lb_coarse_r1.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/synapse/gold_standards/lb_coarse_r2.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/synapse/gold_standards/lb_coarse_r3.csv',
    ], 'fine':[
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/synapse/gold_standards/lb_fine_r1.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/synapse/gold_standards/lb_fine_r2.csv',
        '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/synapse/gold_standards/lb_fine_r3.csv',
    ]}
    
    def __init__(self):
        self.gold = dict()
        for grain, paths in self.gs_paths.items():
            dfs = list()
            for path in paths:
                df = pd.read_csv(path).set_index(['dataset.name', 'sample.id', 'cell.type'])
                df.columns = ['measured']
                dfs.append(df)
            self.gold[grain] = pd.concat(dfs)
        
    def plot(self,plotDat, tag=None, log=True, by='cell_type', data_set=None, title=None, alpha=.4):
        pDat = plotDat.copy()
        gcorr = pearsonr(pDat.measured, pDat.prediction)[0]
        corrs = pDat.groupby(pDat[by]).apply(lambda x: pearsonr(x.measured, x.prediction)[0])
        pDat['corr'] = corrs[pDat[by]].values
        by_str = '{}_pearson'.format(by)
        pDat[by_str] = pDat.apply(lambda x: '{} {:.2f}'.format(x[by], corrs[x[by]]), axis=1)
        if data_set:
            pDat = pDat.loc[pDat['dataset_name']==data_set]
        pl = (pn.ggplot(pn.aes('measured', 'prediction', color=by_str), pDat)
             + pn.geom_point(alpha=alpha)
             + pn.stat_smooth(mapping=pn.aes('measured', 'prediction', color=by_str),
                              method='lm', geom='line', alpha=0.5, se=False, inherit_aes=False))
        if len(pDat['sample'].unique()) < 10:
            pl = pl + pn.aes(shape='sample')
        else:
            pl = pl + pn.aes(shape='dataset_name')
        if log is True:
            pl = pl + pn.scale_x_log10() + pn.scale_y_log10()
        if title is not None:
            pl = pl + pn.ggtitle(title)
        elif tag is not None:
            pl = pl + pn.ggtitle('{} pearson={}'.format(tag, gcorr))
        else:
            pl = pl + pn.ggtitle('pearson={}'.format(gcorr))
        return pl

    def make_plot_data(self, silver, grain='coarse'):
        dat = pd.concat([silver.dropna(), self.gold[grain].dropna()], axis=1, join='inner')
        new_index = [re.sub('\.', '_', name) for name in dat.index.names]
        dat.index = dat.index.rename(new_index)
        pDat = dat.reset_index()
        pDat['sample'] = pDat.apply(lambda x: '{} {}'.format(x.dataset_name, x.sample_id), axis=1)
        pDat.rename(columns={'dataset.name': 'dataset_name'}, inplace=True)
        pDat['ds_cell_type'] = pDat.cell_type.str.cat(pDat.dataset_name, sep=' ')
        return pDat
    
    def get_silver(self, tag, grain='coarse'):
        silver_path = '../docker/docker_outs/out_DREAM_' + grain + '_'+tag+'/predictions.csv'
        silver = pd.read_csv(silver_path)
        silver = silver.set_index(['dataset.name', 'sample.id', 'cell.type'])
        return silver

    def show(self, tag, log=True, by='cell_type', grain='coarse', data_set=None, title=None, alpha=.4):
        silver = self.get_silver(tag, grain)
        pDat = self.make_plot_data(silver, grain=grain)
        return self.plot(pDat, tag, log, by, data_set, title=title, alpha=alpha)

    def cor_per_run(self, pattern=r'.*sy', by='ds_cell_type', path='../docker/docker_outs/', grain='coarse'):
        files = os.listdir(path)
        def lines():
            for fl in tqdm(files):
                if re.match(pattern, fl) is None or grain not in fl:
                    continue
                csv = os.path.join(path, fl, 'predictions.csv')
                if not os.path.isfile(csv):
                    continue
                try:
                    silver = pd.read_csv(csv)
                    silver = silver.set_index(['dataset.name', 'sample.id', 'cell.type'])
                    pDat = self.make_plot_data(silver, grain=grain)
                except:
                    continue
                sm = lambda x: spearmanr(x.measured, x.prediction, nan_policy='omit')[0]
                pm = lambda x: pearsonr(x.measured, x.prediction)[0]
                scorrs = pDat.groupby(pDat[by]).apply(sm)
                pcorrs = pDat.groupby(pDat[by]).apply(pm)
                yield {
                    'tag':fl.split(f'DREAM_{grain}_')[-1],
                    'spear_mean':np.mean(scorrs),
                    'spear_median':np.median(scorrs),
                    'pear_mean':np.mean(pcorrs),
                    'pear_median':np.median(pcorrs),
                }
        return pd.DataFrame(lines()).set_index('tag').sort_values(by='pear_mean', ascending=False)
    
    def loss_hist(self, tag, grain='coarse', smooth=100, base='../docker/docker_outs/out_DREAM'):
        path = base + '_' + grain + '_' + tag
        def load_one_hist(fl, path=path, smooth=smooth):
            if re.match(r'approx_.*', fl) is None:
                return
            pkl = os.path.join(path, fl)
            with open(pkl, 'rb') as buff:
                hist = pickle.load(buff).hist
            if any(hist<=0):
                shist = gaussian_filter1d(hist, sigma=smooth)
            else:
                shist = np.exp(gaussian_filter1d(np.log(hist), sigma=smooth))
            m = re.search('approx_(.+)_(.+).pkl', fl)
            return pd.DataFrame({'dataset': m.group(1),
                               'sample': m.group(2),
                               'id': m.group(1) + ' ' + m.group(2),
                               'iteration': range(len(hist)),
                               'hist': hist,
                               'smooth hist': shist})
        def hists(path, smooth=smooth):
            files = os.listdir(path)
            with mp.Pool() as pool:
                for sdf in tqdm(pool.imap(load_one_hist, files),
                                total=len(files), desc='loading loss history'):
                    yield sdf
        return pd.concat(hists(path))


