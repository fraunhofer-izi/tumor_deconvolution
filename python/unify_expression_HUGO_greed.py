import os
import sys
import argparse
import pickle
import re
import multiprocessing as mp
from multiprocessing.reduction import ForkingPickler, AbstractReducer
import rpy2
from tqdm.auto import tqdm
import pandas as pd
import numpy as np
from scipy.special import logsumexp
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
make_names = robjects.r['make.names']
os.environ['PYTHONHASHSEED'] = '0' # no salt in hash function


# +
class ForkingPickler4(ForkingPickler):
    def __init__(self, *args):
        if len(args) > 1:
            args[1] = 2
        else:
            args.append(2)
        super().__init__(*args)

    @classmethod
    def dumps(cls, obj, protocol=4):
        return ForkingPickler.dumps(obj, protocol)

def dumper(obj, file, protocol=4):
    ForkingPickler4(file, protocol).dump(obj)

class Pickle4Reducer(AbstractReducer):
    ForkingPickler = ForkingPickler4
    register = ForkingPickler4.register
    dump = dumper

ctx = mp.get_context()
ctx.reducer = Pickle4Reducer()
# -

if os.getcwd().endswith('python'):
    os.chdir('..')

if __name__ == '__main__' and not sys.argv[-1].endswith('json'):
    desc = 'Obtimizes feature selection based on availability and mapping quality.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--start', help=f'Start with line sperated features in passed file.',
                        type=str, default=None, metavar='feature_file.txt')
    parser.add_argument('--score', help='Group features by score instead of sample signature.',
                        action='store_true')
    parser.add_argument('--onlySeq', help='Use only sequencing data.',
                        action='store_true')
    parser.add_argument('--remove', help='Also allow to remove features during optimization.',
                        action='store_true')
    parser.add_argument('--tag', help='Output file name tag.',
                        type=str, default='opti', metavar='tag')
    parser.add_argument('--scale', help='A factor to dowscale the score vector befor applying the softmax.',
                        type=np.float, default=10, metavar='s')
    parser.add_argument('--patience', help='Number of iterations to wait for score to increase before stopping.',
                        type=int, default=100, metavar='integer')
    parser.add_argument('--seqw', help='Weight factor for RNA-Seq data relative to array data.',
                        type=np.float, default=2, metavar='float')
    parser.add_argument('grain', help='"coarse" or "fine" grain',
                        metavar='grain', type=str, nargs='?', default='coarse')
    args = parser.parse_args()
elif sys.argv[-1].endswith('json'):
    # We are probably in a jupyter notebook
    class A:
        start = '/<censored_path>/dominik.otto/tumor-deconvolution-dream-challenge/feature_whseq_selection_fine.txt'
        grain = 'fine'
        tag = 'interactive'
        remove = False
        score = False
        scale = 10
        patience = 10
        seqw = 5
        onlySeq = True
    args = A()

with_array = True
rescale_arrays = True
rescale_array_to_n_reads = 5e7
quality_threshold = 5
quality_default = 7
exlcuded_cohorts = [
    'GSE79706', # most monocytes seem like macrophages
    'GSE41914', # only cohort wrongly classified as naive B-cells while labled as NK-cells
    'GSE116672', # viscRNA-Seq that are very inconsistent
    'GSE66360', # by far most wrongly classified as CD4 T-cells while labled as endothelia cells
    'GSE31312' # by far most wrongly classified as endothelial cells while labled as B-cells
]
softmax_downscale = args.scale
important_genes_table = "tabula-SupplementaryTables.csv"
importance = 2 # weight factor for important genes
workPath = "../data"
workPath = os.path.abspath(workPath)
inRDS = os.path.join(workPath, "rnaseq-R-HUGO-GSE")
inRDSarray = os.path.join(workPath, "rnaseq-R-HUGO-GSE-arrays")
outName = "tissue_tsv_greed"
if with_array:
    outName += "_array"
outTSV = os.path.join(workPath, outName)
HUGOfile = os.path.join(workPath, "HUGO_mpa.RDS")
fl_name = os.path.join(inRDS, "info.RDS")
rnaseq_fl_name = os.path.join(inRDS, "info.RDS")
array_fl_name = os.path.join(inRDSarray, "info.RDS")
annoFile = os.path.join(workPath, "sample_data.tsv")
grain = args.grain
report_out = os.path.join(workPath, f'feature_{args.tag}_report_{grain}.pkl')
selection_out = os.path.join(workPath, f'feature_{args.tag}_selection_{grain}.txt')
if args.start:
    report_out = os.path.join(workPath, f'feature_{args.tag}_report_{grain}_start.pkl')
    selection_out = os.path.join(workPath, f'feature_{args.tag}_selection_{grain}_start.txt')

try:
    os.mkdir(outTSV)
except FileExistsError:
    pass

print(f'grain {args.grain} start {args.start} score {args.score}')

robjects.r('''
    info_combine = function(A, B) {
        if (!is.list(A) || !is.list(B)) return(c(A, B))
        names(A) = make.names(names(A))
        names(B) = make.names(names(B))
        info = list()
        for (k in union(names(A), names(B))) {
            if (!k %in% names(A)) {
                info[[k]] = B[[k]]
            } else if (!k %in% names(B)) {
                info[[k]] = A[[k]]
            } else {
                info[[k]] = info_combine(A[[k]], B[[k]])
            }
        }
        return(info)
    }
''')
info_combine = robjects.r['info_combine']
readRDS = robjects.r['readRDS']

HUGO_map = readRDS(HUGOfile)
HUGO_map = pandas2ri.rpy2py_dataframe(HUGO_map)
imp_genes = pd.read_csv(important_genes_table, comment="#", low_memory=False)
anno = pd.read_csv(annoFile, sep='\t', low_memory=False)

selection = set()
if args.start:
    with open(args.start, 'r') as f:
        selection = set(line.rstrip('\n') for line in f.readlines())

features = set(HUGO_map['HUGO']) - {'none'}

print('Loading and combining Rlang mapping metadata ...', file=sys.stderr)
r_infos = info_combine(readRDS(array_fl_name), readRDS(rnaseq_fl_name))
print('done.', file=sys.stderr)

# getting important genes of The Cancer Immunome Atlas (https://tcia.at/)
imp_genes['Metagene'] = imp_genes['Metagene'].str.replace('ClQA', 'C1QA')
imp_genes['Metagene'] = imp_genes['Metagene'].str.replace('ClQB', 'C1QB')
imp_genes['Metagene'] = imp_genes['Metagene'].str.capitalize()

HUGO_map['cap_ref'] = HUGO_map['ref'].str.capitalize()
refs = set(HUGO_map['cap_ref'])
missing = {gene for gene in imp_genes['Metagene'] if gene not in refs}
if missing:
    print(f'The following {len(missing)} genes of "The Cancer Immunome Atlas" '
          f'were not found in the HUGO genes: {list(missing)}', file=sys.stderr)
important_genes = refs - missing
important_genes = set(HUGO_map.set_index('cap_ref').loc[important_genes, 'HUGO']) - {'none'}


def normalize_array(expr):
    expr -= expr.min()
    if rescale_arrays:
        expr *= rescale_array_to_n_reads / expr.sum()
    return expr


def filtered_expr(set_name, grain, cell_type, features):
    set_name = set_name.replace(".p", "-p")
    is_array = re.compile(r'^[^-]*.p[0-9]*$').match(set_name)
    if is_array:
        inDir = inRDSarray
    else:
        inDir = inRDS
    fl_name = os.path.join(inDir, '-'.join([set_name, grain, make_names(cell_type)[0] + '.RDS']))
    expr = pandas2ri.rpy2py_dataframe(readRDS(fl_name))
    expr = expr.reindex(features)
    if is_array:
        expr = normalize_array(expr)
    return expr


def filtered_samples_expr(set_name, grain, cell_type, features):
    expr = filtered_expr(set_name, grain, cell_type, features)
    local_gse = re.sub(r'-.*', '', set_name)
    local_anno = anno.loc[(anno['gse'] == local_gse) & ~anno['aggregated_quality'].isna(), :]
    quality = local_anno.drop_duplicates('id').set_index('id').reindex(set(expr.columns))
    quality = quality.loc[expr.columns]['aggregated_quality'].fillna(quality_default)
    return expr.loc[:, quality > quality_threshold]


def get_availability_frame(info, seq_weight=args.seqw, onlySeq=args.onlySeq):
    grain, cell_type, set_name, r_df = info
    for gse in exlcuded_cohorts:
        if gse in set_name:
            return None
    is_array = re.compile(r'^[^-]*.p[0-9]*$').match(set_name)
    if is_array and onlySeq:
        return None
    info_df = pandas2ri.rpy2py_dataframe(r_df).set_index('HUGO').reindex(features)
    total = info_df['total'][0]
    replacements = dict()
    for col in columns:
        if col in na2total:
            replacements[col] = total
        else:
            replacements[col] = 0
    info_df.fillna(replacements, inplace=True)
    expr = filtered_samples_expr(set_name, grain, cell_type, features)
    good_feature = info_df['hits'] - info_df['min_missing'] >= 1
    good_feature &= info_df['max_ambiguity'] == 1
    good_sample = 1 if is_array else seq_weight
    return (~expr.isna()).multiply(good_feature + weight, axis=0) * good_sample


infos = dict()
columns = {'hits', 'max_ambiguity', 'min_missing', 'NAs', 'total'}
na2total = {'min_missing', 'NAs', 'total'}
uninteresting = {'too.unspecific', 'NA.', 'others'}
weight = pd.Series(1, index=features)
weight[important_genes] = 2
for g, ct_list in r_infos.items():
    if g == grain:
        break
avail = dict()
for cell_type, ds_list in tqdm(ct_list.items(), total=len(ct_list), desc=grain):
    if cell_type in uninteresting:
        continue
    tasks = [(grain, cell_type) + item for item in ds_list.items()]
    avails = list()
    with mp.Pool() as pool:
        for expr in tqdm(
            pool.imap_unordered(get_availability_frame, tasks),
            total=len(tasks),
            desc=cell_type,
            disable=True
        ):
            if expr is not None:
                avails.append(expr)
    avail[cell_type] = pd.concat(avails, axis=1)

for ct, df in avail.items():
    print(f'{ct}: {df.shape[1]}', file=sys.stderr)


def universal_get_score(selection, avail=avail):
    scores = dict()
    all_samps = set()
    for ct, expr in avail.items():
        local_expr = expr.loc[selection, :]
        ind = (local_expr > 0).all()
        scores[ct] = local_expr.loc[:, ind].values.sum()
        all_samps.update(ind[ind].index)
    scores['total'] = -logsumexp(-np.array(list(scores.values()))/softmax_downscale)
    scores['signature'] = hash(frozenset(all_samps))
    return scores


iteration = 0
bad_iterations = 0
report = list()
feature_avail = dict()
if args.score:
    group_stat = 'total'
else:
    group_stat = 'signature'
if selection:
    current = universal_get_score(selection)
    current_score = current['total']
    current_signature = current['signature']
else:
    current_score = None
    current_signature = None

for f in tqdm(features, desc='prepare features'):
    feature_avail[f] = dict()
    for ct, df in avail.items():
        feature_avail[f][ct] = df.loc[f, :]


def get_sample_wise(expr, selection):
    local_expr = expr.loc[selection, :]
    ind = (local_expr > 0).all()
    return local_expr.loc[:, ind].sum()


def get_score(info):
    feature, selection, avail, block = info
    if feature in selection:
        local_sel = selection.copy()
        local_sel.remove(feature)
        scores = universal_get_score(local_sel, avail)
        scores['remove'] = True
    else:
        scores = dict()
        all_samps = set()
        for ct, vec in block.items():
            new_vec = avail[ct]
            new_vec = new_vec[new_vec > 0]
            samples = set(vec.index).intersection(new_vec.index)
            all_samps.update(samples)
            scores[ct] = vec.loc[samples].values.sum() + new_vec.loc[samples].values.sum()
        scores['total'] = -logsumexp(-np.array(list(scores.values()))/softmax_downscale)
        scores['signature'] = hash(frozenset(all_samps))
        scores['remove'] = False
    scores['feature'] = feature
    return feature, scores


def list_to_text(a_list, n_max=10):
    if len(a_list) == 0:
        return '-'
    elif len(a_list) == 1:
        return str(a_list[0])
    elif len(a_list) <= n_max:
        return str(a_list)
    else:
        text = '['
        for i in range(n_max):
            text += str(a_list[i]) + ', '
        text += '... ]'
        return text


if selection:
    print(f'Starting with {len(selection)} genes and score {current_score}.',
          file=sys.stderr)
else:
    print('Starting with 0 genes.', file=sys.stderr)
while True:
    new_scores = dict()
    block = dict()
    for ct, df in tqdm(avail.items(), desc='prepare', disable=True):
        block[ct] = get_sample_wise(df, selection)
    sel_avail = dict()
    for ct, df in avail.items():
        sel_avail[ct] = avail[ct].loc[selection, :]
    def task_iter(features):
        for f in features:
            if f in selection:
                if args.remove:
                    yield f, selection, sel_avail, None
            else:
                yield f, selection, feature_avail[f], block
    with mp.Pool() as pool:
        if args.remove:
            ntotal = len(features)
        else:
            ntotal = len(features) - len(selection)
        feature_iter = tqdm(
            pool.imap_unordered(get_score, task_iter(features), 100),
            total=ntotal,
            desc=f'itertaion {iteration}'
        )
        for feature, score in feature_iter:
            new_scores[feature] = score
    result_df = pd.DataFrame(new_scores).T
    result_df = result_df.loc[(result_df['signature'] != current_signature) | (result_df['remove'] == False), :]
    group = result_df.loc[result_df['total'].astype(np.float).idxmax(), group_stat]
    result = result_df.loc[result_df[group_stat] == group, :]
    next_selection = selection.copy()
    next_features = set(result['feature'])
    if result['remove'][0]:
        next_selection -= next_features
        diff = -len(result.index)
    else:
        next_selection.update(next_features)
        diff = len(result.index)
    total = universal_get_score(next_selection)
    total['feature'] = 'total'
    result = result.append(pd.DataFrame(total, index=[1]))
    next_score = total['total']
    nfeatures = len(next_selection)
    result.loc[:, 'nfeatures'] = nfeatures
    result.loc[:, 'iteration'] = iteration
    report.append(result)
    genes = list_to_text(list(next_features))
    prob = nfeatures / len(features)
    print(f'genes: {nfeatures} ({prob:.2%}), diff: {diff}, score: {next_score}, genes: {genes}',
          file=sys.stderr)
    report_df = pd.concat(report).reset_index()
    with open(report_out, 'wb') as f:
        pickle.dump(report_df, f)
    selection = next_selection
    if current_score and (next_score <=  current_score):
        bad_iterations += 1
        print(f'Bad itertaion {bad_iterations}.', file=sys.stderr)
    else:
        bad_iterations = 0
        current_score = next_score
        with open(selection_out, 'w') as f:
            for feature in selection:
                f.write(feature + '\n')
    if bad_iterations >= args.patience:
        print("Converged.", file=sys.stderr)
        break
    iteration += 1
    if len(selection) == len(features):
        print('All selected.', file=sys.stderr)
        break


