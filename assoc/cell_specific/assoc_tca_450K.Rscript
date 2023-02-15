#!/usr/bin/env Rscript --vanilla --slave
##################################################################################################
# R script for running linear models for bulk tissue using limma (https://bioconductor.org/packages/release/bioc/html/limma.html)
# and cell-specific analysis using TCA (https://github.com/cozygene/TCA)
# day: 10 February 2023
# author: Jose Jaime Martinez-Magana
####################################################################################

# This script uses three inputs the qcdata following this github https://github.com/martinezjaime/ewas_saliva_ses/blob/main/qc_data/probe_and_sample_quality_control.ipynb,
# a list of samples to be included in the analysis as csv file
# the path for the results of the epistructure following this github: https://github.com/martinezjaime/ewas_saliva_ses/blob/main/epigenetic_ancestry/epigenetic_ancestry.ipynb

####################################################################################
# set parameters
# this script uses the library optparse to add arguments to the script
# adding arguments to the script
library(optparse) 
option_list=list(
    make_option(c("--file"),
                type="character",
                default=NULL,
                help="path to *rds object with qc information", metavar="character"),
    make_option(c("--out"),
                type="character",
                default="out.rds",
                help="output file name for rds [default= %default]", metavar="character"),
    make_option(c("--outname"),
                type="character",
                default=NULL,
                help="name of the rds file after performing this analysis", metavar="character"),
    make_option(c("--samplelist"),
                type="character",
                default=NULL,
                help="path to a *csv file with a list of samples to subset the analysis. This script uses a file with SampleID with Array_Sentrix structures and a header with SampleID", metavar="character"),
    make_option(c("--pcafile"),
                type="character",
                default=NULL,
                help="path to pca output from epistructure-glint for epigenetic ancestry", metavar="character"),
    make_option(c("--npcs"),
                type="integer",
                default=2,
                help="Add the number of ancestry PCs to adjust the models, default = 2", metavar="character"),
    make_option(c("--pheno"),
                type="character",
                default=NULL,
                help="name of the column identifying phenotype to be tested in the regression model", metavar="character"),
    make_option(c("--covar"),
                type="character",
                default=NULL,
                help="list of covariates to be included in the analysis from the phenofile separeted by ','", metavar="character"),
    make_option(c("--sva"),
                type="character",
                default=NULL,
                help="add TRUE if you want the model to be adjusted by subrrogate variables using SVA", metavar="character"),
    make_option(c("--nsva"),
                type="integer",
                default=2,
                help="Add the number of sva to adjust the models, default = 2", metavar="character"),
    make_option(c("--mprop"),
                type="character",
                default=NULL,
                help="add TRUE if you want the model to be adjusted by proportions of DNA methylation", metavar="character")
);

opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

####################################################################################
# load libraries
library(preprocessCore)
library(lumi)
library(TCA)
library(minfi)
library(pracma)
library(matrixStats)
library(limma)
library(sva)
library(SmartSVA)
library(magrittr)
library(plyr)
library(dplyr)


# loading rds
paste0("Start analysis of data:",Sys.time(),"---","###Analysis path[",opt$file,"]###")
file=readRDS(opt$file)
# setting output file
outfile=paste0(opt$out,opt$outname,sep="")
paste0("Output of cell-type specific statistical models will be saved to:",Sys.time(),"---","###Output path[",outfile,"]###")

if(!is.null(opt$pcafile)){
    # loading epigenetic ancestry pcs file
    paste0("Loading epigenetic ancestry:",Sys.time())
    paste0("Warning: this script will use any genetic ancestry estimated by the user")
    pcs=read.table(opt$pcafile)
    # add colnames to pcs
    colnames(pcs)=c("SampleID",paste(rep('ances_pc',10),rep(1:10),sep=''))
    # add rownames
    rownames(pcs)=pcs$SampleID
    # remove innecesary column
    pcs$SampleID = NULL
    # subsetting pcs 
    paste0("Warning: this script will use the number of epigenetic ancestry principal components specify, default = 2")
    # subsetting the number of pcs to user identified
    pcs=pcs[,c(1:opt$npcs)]
}

# loading phenotype file
paste0("Loading pheno file:",Sys.time())
# add pheno to be tested
pheno=opt$pheno
# add covariates to be tested
covar=unlist(strsplit(opt$covar,','))
# get covars from pheno file
covars=c(pheno,covar)
covars=c(pheno,covar,"SampleID")
# select covars from phenofile
pheno_data=file$pheno[,colnames(file$pheno) %in% covars]
# add SampleID as rownames
rownames(pheno_data)=pheno_data$SampleID
# remove SampleID column
pheno_data$SampleID=NULL

# merge pheno_data with epigenetic ancestry pcs
paste0("Merging pheno file with epigenetic ancestry:",Sys.time())
pheno_data=cbind(pheno_data,pcs)

# get cells
paste0("Adding cells to pheno file:",Sys.time())
cells=file$cell_epidish
# transform cells to data frame
cells2=as.data.frame(cells)
# remove unneccesary object
rm(cells)
# removing cells with cero numbers in all samples
cells3=cells2[,colSums(cells2) > 0.00]
paste0("Warning: This script merges CD4T + CD8T cells to Limphocytes imputed with epiDISH from saliva, please revised the cell distribution in your data")
# merging CD4T and CD8T cells
cells3_m=cells3[,"CD4T", drop = F] + cells3[,"CD8T", drop = F]
colnames(cells3_m)="Lympho"
cells3_m2=cbind(cells3_m, cells3[,c("Epi","Fib","B","NK","Mono","Neutro")])
# add warning to say the user the cells that will be used in the analysis
paste0("Using the following cells in the regression models:",colnames(cells3_m2))
# merge pheno_data with cells
pheno_data=cbind(pheno_data,cells3_m2)

if(opt$mprop){
    paste0("Adding proportion of methylated and unmethylated pheno file:",Sys.time())
    # get methylated and unmethylated proportions as covariates
    umm=file$covariates_umm
    # merge pheno_data with methylated and unmethylated proportions
    pheno_data=cbind(pheno_data,umm)
} else {
    paste0("Methylation proportions TURNED OFF not adjusting for this covariates",Sys.time())
    pheno_data=pheno_data
}

######################################################################################
# subset samples if sample list is added
# subset phenofile
# here we used is.null, to ask if there was a sample list added to the script
if(!is.null(opt$samplelist)){
    paste0("Sample filter added, phenofile will be filtered based in:",opt$samplelist,"---",Sys.time())
    sample_f=read.csv(opt$samplelist)
    # subsetting pheno file
    pheno_data=pheno_data[rownames(pheno_data) %in% sample_f$SampleID,]
} else {
    pheno_data=pheno_data
}

# subset beta matrix
if(!is.null(opt$samplelist)){
    paste0("Sample filter added, beta matrix will be filtered based in:",opt$samplelist,"---",Sys.time())
    sample_f=read.csv(opt$samplelist)
    bvalues=file$rgsetraw_bmiq
    bvalues=bvalues[,colnames(bvalues) %in% sample_f$SampleID]
} else {
    bvalues=file$rgsetraw_bmiq
}

# sample sanity check
all_samples=is.null(pheno_data[!colnames(bvalues) %in% rownames(pheno_data),])
if(all_samples == TRUE){
    paste0("All samples in the phenofile found in the beta value matrix")
    bvalues=bvalues
} else {
    missing_samples=rownames(pheno_data[!rownames(pheno_data) %in% colnames(bvalues),])
    pheno_data=pheno_data[!rownames(pheno_data) %in% missing_samples,]
    bvalues=bvalues[, colnames(bvalues) %in% rownames(pheno_data)]
    paste0("The following samples were not found in the beta values but found in the pheno file: ", missing_samples)
}

# ordering pheno file based on beta
paste0("Warning: this script orders the pheno file based on the beta matrix:", Sys.time())
pheno_data=pheno_data[match(rownames(pheno_data),colnames(bvalues)),]

######################################################################################
# transform beta matrix to m-values using lumi package for statistical testing
mvalues=beta2m(bvalues)

######################################################################################
# estimate subrrogate variables using SmartSVA
if(opt$sva){
    pheno_sva=opt$pheno
    # create model matrix with phenotype
    mod=model.matrix(~pheno_data[,pheno_sva], data=pheno_data)
    # remove the effect of the phenotype in the mvalue matrix
    # this matrix has to be transpose two times
    mvalues_r=t(resid(lm(t(mvalues) ~ pheno_data[,pheno_sva], data=pheno_data)))
    # create null model
    # this command is needed for SVA estimating using sva package, mod0=model.matrix(~1,data=pheno_data)
    # estimating number of subrrogate variables
    paste0("Computing subrrogate variables with SmartSVA:", Sys.time())
    n_sv=EstDimRMT(mvalues_r, FALSE)$dim + 1
    # this command is needed for SVA estimating using sva package, n_sv=num.sv(mvalues,mod,method="leek")
    svobj=smartsva.cpp(mvalues,mod,mod0=NULL,n.sv=n_sv)
    # adding subrrogate variable to phenofile
    svobj_df=data.frame(svobj$sv)
    # adding sampleIDs to rownames
    rownames(svobj_df)=rownames(pheno_data)
    # adding colnames
    colnames(svobj_df)=c(paste0(rep('sva',n_sv),rep(1:n_sv)))
    paste0("Warning: this analysis will be adjusted by the specify number subrrogate variables, default = 2")
    paste0("Warning: adjust this number if needed")
    # subsetting the number of sva to user identified
    svobj_df=svobj_df[,c(1:opt$nsva)]
    # merging with pheno file
    pheno_data=cbind(pheno_data, svobj_df)
    # ordering pheno file based on beta values
    paste0("Warning: this script orders again the pheno file based on the beta matrix:", Sys.time())
    pheno_data=pheno_data[match(rownames(pheno_data),colnames(bvalues)),]
    paste0("End of computation of subrrogate variables with SmartSVA:", Sys.time())
} else {
    paste0("The sva method is turned off, NO adjustment of the models for sva", Sys.time())
    pheno_data=pheno_data
}

######################################################################################
# this function will run cell-specific analysis using TCA Tensor Composition Analysis
# follow this github for more information: https://github.com/cozygene/TCA
# TCA function requieres four inputs
# X = beta value matrices
# W = cell-type proportions all summing to one, very important
# C1 = covariates to be tested, including your phenotype
# C2 = other confounding covariates into your models, could be PCS, slide, array, plate, etc
# extract cells from phenofile
# this is the W argument for the TCA function
cells_tca=pheno_data[,c("Epi","Fib","B","NK","Mono","Neutro","Lympho")]
paste0("Testing the following cells in TCA: ", colnames(cells_tca))
# this is the C2 argument for the TCA function
conf_tca=svobj_df
# adding slide and array as confounding covariates into TCA
conf_tca$slide=strsplit(rownames(pheno_data), "_") %>% sapply(extract2, 1)
conf_tca$array=strsplit(rownames(pheno_data), "_") %>% sapply(extract2, 2)
# we are going to transform the slide in numeric to fit the linear models
conf_tca$slide=as.numeric(conf_tca$slide)
# we are going to transform the array in numeric to fit the linear models
conf_tca$array=as.numeric(as.factor(conf_tca$array))
paste0("Using the following counfunding covariates in TCA: ", colnames(conf_tca))
paste0("Warning: We are going to adjust the models with SVAs, estimated using SmartSVA")
# subsetting to testing covariates
paste0("Adjusting TCA by the following confounders: ", colnames(conf_tca))
# this is the C1 argument for the TCA function
filter_cov=c(colnames(cells_tca),colnames(conf_tca))
other_covars=pheno_data[!colnames(pheno_data) %in% filter_cov]

# creting an empty list to store output
out=list()
# running tca
tca=tca(X=bvalues,
        W=cells_tca,
        C1=other_covars,
        C2=conf_tca)
# add data to output
out$tca=tca

######################################################################################
# save output
saveRDS(file=outfile,out)
