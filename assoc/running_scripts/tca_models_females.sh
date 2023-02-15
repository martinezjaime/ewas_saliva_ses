#!/bin/bash
#SBATCH --job-name=running_tca_models_ewas_saliva_ses
#SBATCH --out="slurm-%j.out"
#SBATCH --time=04:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --partition=day
####################################################################################
# script for running the TCA models for saliva SES in children
# day: 15 february 2023
# analyzer: Jose Jaime Martinez-Magana
# cluster: Grace - HPC Yale
####################################################################################
# This script uses the assoc_tca_450K.Rscript for linear models that could be found in the 
# following github: https://github.com/martinezjaime/ewas_saliva_ses/tree/main/assoc/cell_specific
####################################################################################
# load miniconda
module load miniconda
# activate environment
conda activate ewas_saliva
# path to Rscript
# this path should be change based on datapaths
Rscript="/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/script/assoc/cell_specific/assoc_tca_450K.Rscript"

# running script for female unrelated individuals
Rscript $Rscript --file=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/databases/qced/qced_data_v02062023.rds \
--out=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/cell_specific/ \
--outname=tca_unrelated_female_v15022023.rds \
--samplelist=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/samplelist/sample_subset_unrelated_females_v02142023.csv \
--pcafile=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/pcs/glint_epi_children_v08February2023.epistructure.pcs.txt \
--pheno=SES \
--covar=age \
--sva=TRUE \
--mprop=FALSE \
--nsva=2 \
--npcs=2
