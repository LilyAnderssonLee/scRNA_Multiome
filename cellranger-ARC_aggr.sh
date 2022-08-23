#!/bin/bash -l

#SBATCH -A snic2021-22-916
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 96:00:00
#SBATCH -J cellranger-ARC-agg

module load bioinfo-tools
module load cellranger-ARC/1.0.0
module load cellranger-ARC-data/2020-A

transcrip=$CELLRANGER_ARC_DATA/refdata-cellranger-arc-mm10-2020-A
csv_file=/proj/snic2021-23-715/private/Lili/analysis/NGI_proj/multiome/Ana/library_csv/CP_aggr.csv

cellranger-arc aggr --id=CP_aggr \
	--reference=$transcrip \
	--csv=$csv_file \
	--normalize=depth \
	--localcores=10 \
	--localmem=60
