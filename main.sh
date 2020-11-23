#! /bin/bash -

echo "Downloading..."
./download-rnaseq-studies.sh
echo "Downloading annotations..."
#./download-rnaseq-studies-soft.sh
./download-studies-soft.sh
echo "Unpacking..."
./unpack.sh
echo "Extracting meta info..."
./make-soft-dict.sh

echo "Resolving all studies ..."
source /etc/profile.d/modules.sh
module load R/3.5.2-0
./resolve_all.sh

echo "Making extraction statistics and test data set ..."
R/extraction_stats.R

echo "Making annotation table ..."
R/annotation_table.R

# Some data set tagging happened here with https://ribogit.izi.fraunhofer.de/Dominik/geotag

echo "Aggregating tagging information..."
source souce ~dominik.otto/enterEnv.sh # conda environment with python 3.6
python3 python/aggregate_tagging_data.py

echo "Downloading and extracting gene info..."
./download_gene_info.sh

echo "Preparing HUGO gene symbol map..."
R/make_HUGO_map.R

echo "Map to HUGO gene symbols ..."
R/map2HUGO.R

echo "Downloading and mapping all microarrays..."
R/get_arrays.R

echo "Selecting features based on availability..."
sbatch slurm_job_feature_selection.sh --tag wsy coarse # Sub_Challenge_1
sbatch slurm_job_feature_selection.sh --tag wsy fine # Sub_Challenge_2

echo "Collecting counts per cell type ..."
R/unify_expression_HUGO_by_selection.R &> >(tee unify_expression_wsy.err)

echo "Making characterization..."
source ~dominik.otto/enterEnv.sh # conda environment with python 3.6
outDir="$(readlink -f ../data)"
tsvinsel="$outDir/tissue_tsv_wsy_array/"
sbatch slurm_job_char.sh --center --simp coarse "$tsvinsel" \
    "$outDir/genchar_coarse_unclean_center_simp_wsy.pkl" # Sub_Challenge_1
sbatch slurm_job_char.sh --center --simp fine "$tsvinsel" \
    "$outDir/genchar_fine_unclean_center_simp_wsy.pkl" # Sub_Challenge_2

echo "Now run python/make_cell_type_characterisation_tagged_wsy_fine.ipynb" # Sub_Challenge_1
echo "and make_cell_type_characterisation_tagged_wsy_coarse.ipynb" # Sub_Challenge_2
echo "in jupyter to make the final cell type charachteristics."

echo "Finally you can generate a deconvolution docker app. See docker/Readme for details."
