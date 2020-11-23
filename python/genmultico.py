#! /usr/bin/env python
#SBATCH --job-name="genmultico"
#SBATCH --output=../logs/slurm.genmultico.%j.out
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=230G
#SBATCH --time=10-00:00:00

import sys
import numpy as np
import pickle
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from genchar import *

if __name__ == "__main__":
    grain = sys.argv[1]
    tsv_in = sys.argv[2]
    out_tag = sys.argv[3]
    print('Starting grain {} tag {}...'.format(grain, out_tag))
    expr = Expressions(grain=grain, TSV_dir=tsv_in)
    ch = MultiCoCharacterizer(expr, n_components=10, n_batches=100, batch_size_factor=100)
    chars = ch.affine_trans
    variations = ch.remaining_variation
    backgrounds = ch.background_chars
    chars_file_name = ('../data/genchar_{}_{}_multi.pkl'.format(grain, out_tag))
    print('Saving {}...'.format(chars_file_name))
    chars_data = dict({'chars':chars, 'features':ch.expr.features, 'variations':variations,
                       'multico':True, 'backgrounds':backgrounds})
    with open(chars_file_name, 'wb') as buff:
        pickle.dump(chars_data, buff, protocol=4)
    trans_file_name = ('../data/transformations_{}_{}_multi.pkl'.format(grain, out_tag))
    print('Saving {}...'.format(trans_file_name))
    trans_data = dict({'expr':expr, 'ch':ch})
    with open(trans_file_name, 'wb') as buff:
        pickle.dump(trans_data, buff, protocol=4)
