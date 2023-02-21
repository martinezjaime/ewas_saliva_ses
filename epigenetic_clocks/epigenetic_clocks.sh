#!/bin/bash
#SBATCH --job-name=running_epigenetic_clocks
#SBATCH --out="slurm-%j.out"
#SBATCH --time=04:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --partition=day
####################################################################################
# script for running epigenetic clocks in children
# day: 21 february 2023
# analyzer: Jose Jaime Martinez-Magana
# cluster: Grace - HPC Yale
####################################################################################
# This script uses the epigenetic_clock.Rscript for estimating epigenetic clocks
# following github: https://github.com/martinezjaime/ewas_saliva_ses/tree/main/epigenetic_clocks
####################################################################################
# load miniconda
module load miniconda
# activate environment
conda activate ewas_saliva
# path to Rscript
# this path should be change based on datapaths
Rscript="/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/script/epigenetic_age/epigenetic_clock.Rscript"

# running script for unrelated individuals
Rscript $Rscript --file=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/databases/qced/qced_data_v02062023.rds \
--out=/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/epigenetic_age/ \
--outname=epigenetic_clocks_v15022023.rds \
--covar=age,gender
