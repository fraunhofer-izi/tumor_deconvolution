import os
import glob
import pandas as pd
import numpy as np
import yaml
from tqdm.auto import tqdm

tag_path = '/<censored_path>/dominik.otto/geotag_collect/'
taggers = {'michael.rade', 'christoph.kaempf', 'c.schimmelpfennig', 'dominik.otto', 'alex.scholz', 'markus.kreuz', 'kristin.reiche'}
assigner = 'assigner'

tag_files = dict()
for file in glob.glob(tag_path+'/*.yml'):
    name = os.path.basename(file)[:-4]
    if name in taggers:
        tag_files[name] = file
    elif name == assigner:
        pass
    else:
        print(f'The unknown tagger "{name}" is ignored.')
for tagger in taggers:
    if tagger not in tag_files.keys():
        print(f'No file found for tagger "{tagger}".')

taggs = dict()
for name, file in tqdm(tag_files.items(), desc='loading tag files'):
    with open(file, 'r') as buff:
        taggs[name] = yaml.load(buff, Loader=yaml.SafeLoader)

with open(os.path.join(tag_path, assigner + '.yml'), 'r') as buff:
    assignes =  yaml.load(buff, Loader=yaml.SafeLoader)['tags']['assigned']
all_inices = [tuple(id.split('_'))[:2] for id in list(assignes.keys())]

dfs = list()
qualities = pd.DataFrame(data=all_inices, columns=['gse', 'id']).set_index(['gse', 'id'])
tag = 'quality'
for name, data in tqdm(taggs.items(), desc='making data frames'):
    qualities_dict = data['tags'][tag]
    df_data = [tuple(id.split('_'))[:2]+(q,) for id, q in qualities_dict.items()]
    df = pd.DataFrame(df_data, columns =['gse', 'id', f'{tag}_{name}']).set_index(['gse', 'id'])
    df = df.loc[df.index.unique()]
    dfs.append(df)
    qualities = qualities.join(df, how='outer')

def quality_stats(row):
    mean = np.mean(row)
    min = np.min(row)
    std = np.std(row)
    n_qual = np.sum(~row.isna())
    return pd.Series([mean, std, min, n_qual], index=['mean_quality', 'std_quality', 'min_quality', 'n_qualities'])

tqdm.pandas(desc='quality stats per sample')
one_gsm_quality = qualities.reset_index().drop(columns='gse').groupby('id').mean()
stats = one_gsm_quality.progress_apply(quality_stats, axis=1, result_type='expand')
quali_stats = one_gsm_quality.join(stats)

rseq = pd.read_csv('/homes/olymp/michael.rade/playground/geotag/extraction_stats_MiR.tsv', sep='\t')
arr = pd.read_csv('/homes/olymp/michael.rade/playground/geotag/extraction_stats_array_MiR.tsv', sep='\t')
tt = pd.concat([rseq, arr])

ai = tt.join(quali_stats, on='id', how='right')
#ai.loc[ai['n_qualities'].isna(), 'n_qualities'] = 0
ai['n_qualities'] = ai['n_qualities'].astype(int)
ai.to_csv('/<censored_path>/dominik.otto/geo_tagging.tsv_tmp', sep='\t', index=False)
os.rename('/<censored_path>/dominik.otto/geo_tagging.tsv_tmp',
          '/<censored_path>/dominik.otto/geo_tagging.tsv')
print('Saved.')
