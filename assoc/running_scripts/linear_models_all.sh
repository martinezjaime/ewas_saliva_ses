#!/bin/bash
#SBATCH --job-name=running_linear_models_ewas_saliva_ses
#SBATCH --out="slurm-%j.out"
#SBATCH --time=04:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --partition=day
####################################################################################
# script for running the linear models for saliva SES in children
# day: 15 february 2023
# analyzer: Jose Jaime Martinez-Magana
# cluster: Grace - HPC Yale
####################################################################################
# This script uses the assoc_linear_450K.Rscript for linear models that could be found in the 
# following github: https://github.com/martinezjaime/ewas_saliva_ses/tree/main/assoc/bulk_tissue
####################################################################################
# load miniconda
module load miniconda
# activate environment
conda activate ewas_saliva
# path to Rscript
# this path should be change based on datapaths
Rscript="/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/script/assoc/bulk_tissue/assoc_linear_450K.Rscript"

# running script for unrelated individuals
Rscript $Rscript --file=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/databases/qced/qced_data_v02062023.rds \
--out=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/bulk_tissue/ \
--outname=linear_unrelated_v15022023.rds \
--samplelist=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/samplelist/sample_subset_unrelated_v02142023.csv \
--pcafile=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/pcs/glint_epi_children_v08February2023.epistructure.pcs.txt \
--pheno=SES \
--covar=age,gender \
--sva=TRUE \
--mprop=FALSE \
--nsva=2 \
--npcs=2

# running script for female unrelated individuals
Rscript $Rscript --file=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/databases/qced/qced_data_v02062023.rds \
--out=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/bulk_tissue/ \
--outname=linear_unrelated_female_v15022023.rds \
--samplelist=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/samplelist/sample_subset_unrelated_females_v02142023.csv \
--pcafile=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/pcs/glint_epi_children_v08February2023.epistructure.pcs.txt \
--pheno=SES \
--covar=age \
--sva=TRUE \
--mprop=FALSE \
--nsva=2 \
--npcs=2

# running script for males unrelated individuals
Rscript $Rscript --file=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/databases/qced/qced_data_v02062023.rds \
--out=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/bulk_tissue/ \
--outname=linear_unrelated_male_v15022023.rds \
--samplelist=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/samplelist/sample_subset_unrelated_males_v02142023.csv \
--pcafile=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/pcs/glint_epi_children_v08February2023.epistructure.pcs.txt \
--pheno=SES \
--covar=age \
--sva=TRUE \
--mprop=FALSE \
--nsva=2 \
--npcs=2
