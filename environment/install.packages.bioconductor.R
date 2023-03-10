#!/usr/bin/env Rscript
## install required R bioconductor packages, for ewas_saliva
# install requiered bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c(
			"minfi",
			"IlluminaHumanMethylation450kmanifest",
			"IlluminaHumanMethylationEPICmanifest",
			"wateRmelon",
			"readxl",
			"RPMM",
			"FlowSorted.Blood.450k",
			"FlowSorted.Blood.EPIC",
			"EpiDISH",
			"bacon",
			"sva",
			"AnnotationDbi",
			"org.Hs.eg.db",
			"ENmix",
			"variancePartition",
			"edgeR",
			"IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
			"BiocParallel"),
			update = TRUE, ask = FALSE)
