#! /bin/sh

# This script is uses for testing

set -e
testset="$2"
grain="coarse"
tag="DREAM_${grain}_$1$testset"
echo "====== Tag $tag ====="
version="0.8_$tag"

mkdir -p "/<censored_path>/$USER/docker/in_$tag"
mkdir -p "/<censored_path>/$USER/docker/out_$tag"
if [ "$testset" = "gs" ]; then
    cp input_gold/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
    time="0-12:00:00"
elif [ "$testset" = "lgs" ]; then
    cp input_gold/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
    time="7-00:00:00"
elif [ "$testset" = "sy" ]; then
    cp input_synapse/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
    time="0-12:00:00"
elif [ "$testset" = "lsy" ]; then
    cp input_synapse/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
    time="7-00:00:00"
elif [ "$testset" = "sy2" ]; then
    cp input_synapse2/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
    time="0-12:00:00"
elif [ "$testset" = "lsy2" ]; then
    cp input_synapse2/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
    time="7-00:00:00"
elif [ "$testset" = "sy3" ]; then
    if [ "$grain" = "fine" ]; then
        cp input_synapse3f/* "/<censored_path>/$USER/docker/in_$tag/"
    else
        cp input_synapse3/* "/<censored_path>/$USER/docker/in_$tag/"
    fi
    cores=16
    mem=16G
    time="0-12:00:00"
elif [ "$testset" = "lsy3" ]; then
    if [ "$grain" = "fine" ]; then
        cp input_synapse3f/* "/<censored_path>/$USER/docker/in_$tag/"
    else
        cp input_synapse3/* "/<censored_path>/$USER/docker/in_$tag/"
    fi
    cores=16
    mem=16G
    time="7-00:00:00"
elif [ "$testset" = "em" ]; then
    cp input_R2_emulated/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=16
    mem=16G
else
    cp input/* "/<censored_path>/$USER/docker/in_$tag/"
    cores=5
    mem=5G
    time="2-00:00:00"
fi
cp ../python/deconvolute.py ./data/deconvolute.py
if [ "$grain" = "coarse" ]; then
    #cp ../../data/genchar_coarse_fixed_s2_6.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_newsig7793.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_whseq.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_whseqfr.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_wse.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_wsbi.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_wsy.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_ws300.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_wsys.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_wws500.pkl ./data/chars.pkl
    #cp ../../data/genchar_coarse_tagged_wwss.pkl ./data/chars.pkl
    cp ../../data/genchar_coarse_tagged_wws.pkl ./data/chars.pkl
else
    #cp ../../data/genchar_${grain}_HUGO_array.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_whseq_cpr.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wsar.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wsa.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wsbi.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wsyr.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_ws300.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wsysr.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wws500.pkl ./data/chars.pkl
    #cp ../../data/genchar_fine_tagged_wwssr.pkl ./data/chars.pkl
    cp ../../data/genchar_fine_tagged_wwsr.pkl ./data/chars.pkl
fi
#cp ../../data/genchar_${grain}_HUGO_center_fixed.pkl ./data/chars.pkl
#cp ../../data/genchar_${grain}_fixed_s2r.pkl ./data/chars.pkl
cp data/${grain}_grains data/grains
cp ../python/stickbreaking.py ./data/stickbreaking.py

sudo -H ./docker-build \
    ribogit.izi.fraunhofer.de:4567/dominik/tumor-deconvolution-dream-challenge/$grain:$version Dockerfile

sudo docker login ribogit.izi.fraunhofer.de:4567 -u $USER -p $(cat ~/.ssh/gitlab_push_docker)
sudo docker push ribogit.izi.fraunhofer.de:4567/dominik/tumor-deconvolution-dream-challenge/$grain:$version

srun -c $cores --mem="$mem" --time="$time" --job-name="deconv $tag" --partition=pbatch_all \
    ./slurm_run.sh $grain $tag $version

