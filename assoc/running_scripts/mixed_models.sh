#!/bin/bash
#SBATCH --job-name=running_mixed_glint_models_ewas_saliva_ses
#SBATCH --out="slurm-%j.out"
#SBATCH --time=23:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --partition=day
####################################################################################
# script for running the mixed using glint for saliva SES in children
# day: 15 february 2023
# analyzer: Jose Jaime Martinez-Magana
# cluster: Grace - HPC Yale
####################################################################################
# load environment, only needed if miniconda not loaded
module load miniconda
# activate environment for glint
conda activate glint

# set glint path
glint_path="/gpfs/gibbs/project/montalvo-ortiz/jjm262/programs/glint/glint.py"
# input data path, set the same path of the output for the qcdata_to_glint.Rscript
data_p="/vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/databases/qced/glint"

# running mixed linear models with glint for all samples
python ${glint_path} --datafile ${data_p}/qced_data_v02152023_glint_format.glint \
--ewas --lmm --pheno SES \
--covar age gender Epi Fib B NK Mono Neutro Lympho --norm \
--kinship refactor --k 7 \
--out /vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/mixed_models/mixed_related_all_v15022023

# running mixed linear models with glint for males samples
python ${glint_path} --datafile ${data_p}/qced_data_v02152023_glint_format.glint \
--ewas --lmm --pheno SES \
--covar age gender Epi Fib B NK Mono Neutro Lympho --norm \
--kinship refactor --k 7 \
--keep /vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/samplelist/sample_subset_related_males_v02142023.csv \
--out /vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/mixed_models/mixed_related_males_v15022023

# running mixed linear models with glint for females samples
python ${glint_path} --datafile ${data_p}/qced_data_v02152023_glint_format.glint \
--ewas --lmm --pheno SES \
--covar age gender Epi Fib B NK Mono Neutro Lympho --norm \
--kinship refactor --k 7 \
--keep /vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/samplelist/sample_subset_related_females_v02142023.csv \
--out /vast/palmer/scratch/montalvo-ortiz/jjm262/epigenomics/ewas_saliva_ses/results/assoc/mixed_models/mixed_related_females_v15022023
