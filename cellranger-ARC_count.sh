#!/bin/bash -l

#SBATCH -A snic2021-22-916
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 96:00:00
#SBATCH -J cellranger-ARC-count

module load bioinfo-tools
module load cellranger-ARC/1.0.0
module load cellranger-ARC-data/2020-A

sample_id=$1
lib=$2

transcrip=$CELLRANGER_ARC_DATA/refdata-cellranger-arc-mm10-2020-A

cellranger-arc count --id=$sample_id \
--reference=$transcrip \
--libraries=$lib \
--localcores=10 \
--localmem=60

 
