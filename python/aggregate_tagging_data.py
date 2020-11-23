import os
import glob
import pandas as pd
import numpy as np
import yaml
from tqdm.auto import tqdm

tag_path = '/<censored_path>/dominik.otto/geotag_collect/'
taggers = {'michael.rade', 'christoph.kaempf', 'c.schimmelpfennig', 'dominik.otto', 'alex.scholz', 'markus.kreuz', 'kristin.reiche'}
me = 'dominik.otto'
assigner = 'assigner'
annotation_table = '../data/annotation.tsv'
output = '../data/sample_data.tsv'

annotation_df = pd.read_csv(annotation_table, sep='\t', low_memory=False)

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

with open(os.path.join(tag_path, me + '.yml'), 'r') as buff:
    notes =  yaml.load(buff, Loader=yaml.SafeLoader)['tags']['note']
notes = [tuple(id.split('_'))[:2]+(q,) for id, q in notes.items()]
notes = pd.DataFrame(notes, columns =['gse', 'id', 'note']).set_index(['gse', 'id'])
note_df = annotation_df.join(notes, on=notes.index.names)

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

complete_df = note_df.join(quali_stats, on='id', how='right')

# +
# Fixes

# Michael Rade read study for GSE94820 and suggested to use the samples

ind = (complete_df['gse']=='GSE94820') & (complete_df['note']=='monocytes')
complete_df.loc[ind,'cell.type'] = 'Monocytes'
complete_df.loc[ind,'coarse.cell.type'] = 'monocytic.lineage'
complete_df.loc[ind,'fine.cell.type'] = 'monocytes'

ind = (complete_df['gse']=='GSE94820') & (complete_df['note']=='dendritic cells')
complete_df.loc[ind,'cell.type'] = 'Dendritic cells'
complete_df.loc[ind,'coarse.cell.type'] = 'monocytic.lineage'
complete_df.loc[ind,'fine.cell.type'] = 'myeloid.dendritic.cells'

# some studies were tagged inconsistently

complete_df['aggregated_quality'] = complete_df['min_quality']
ind = (complete_df['std_quality']>3) & (complete_df['val'] == 'GCB') # Germinal center B-cells
complete_df.loc[ind, 'aggregated_quality'] = complete_df.loc[ind, 'quality_dominik.otto']
ind = (complete_df['gse'] == 'GSE73502') | (complete_df['gse'] == 'GSE73765') # difficult studies
ind |= complete_df['gse'] == 'GSE87849' # missing wild types
ind |= complete_df['gse'] == 'GSE94820' # some untypical monocytes
complete_df.loc[ind, 'aggregated_quality'] = complete_df.loc[ind, 'quality_dominik.otto']
# -

complete_df.to_csv(output, sep='\t', index=False)


