#! /bin/sh -
#SBATCH --job-name="docker"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
#SBACTH --partition=pbatch_all

grain=$1
tag=$2
version=$3

sudo docker login ribogit.izi.fraunhofer.de:4567 -u $USER -p $(cat ~/.ssh/gitlab_push_docker)
sudo docker pull ribogit.izi.fraunhofer.de:4567/dominik/tumor-deconvolution-dream-challenge/$grain:$version

echo "====== Tag $tag ====="
time sudo -H ./docker-run-generic \
    ribogit.izi.fraunhofer.de:4567/dominik/tumor-deconvolution-dream-challenge/$grain:$version "$tag"

echo "====== Finished $tag ====="

