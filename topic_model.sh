#!/bin/bash -l
#SBATCH	-A	xxx
#SBATCH	-p	core
#SBATCH	-n	7
#SBATCH	-t	4-00:00:00
#SBATCH -J	topic_modeling

module load bioinfo-tools
module load R/4.0.4
module load R_packages/4.0.4
module load miniconda3/4.5.4
module load python3/3.8.7

k=$1

Rscript topic_model.R $k

