######## RNABC pipeline used for data processing in "Using microarray-based clinical subtyping methods for breast cancer in the 
######## era of high-throughput RNA-sequencing" by Pedersen, CB, Nielsen, FC, Rossing, M, Olsen, LR.
######## Code by Christina B. Pedersen <chrbl@bio.dtu.dk> and Lars R. Olsen <lro@bioinformatics.dtu.dk>
######## Updated May 15, 2018


#### Loading required packages ####
libs <- c('citbcmst', 'ggplot2')
for (lib in libs) {
  if (!require(lib, character.only = T)) {
    install.packages(lib, dependencies = T, quiet = T)
  }
} 

source("http://bioconductor.org/biocLite.R")
libs <- c('preprocessCore', 'sva')
for (lib in libs) {
  if (!require(lib, character.only = T)) {
    suppressWarnings(biocLite(lib))
  }
} 

# Load modified CITBCMST script for nicer plots
source('CITBCMST_classifier.R')

#### Loading data ####
# Import full and core CIT data
load('CITBCMST.Rdata')

# Import Bordet data (the data we wish to classify)
load('Bordet.Rdata')


#### Actual data processing ####
### Removal of AFFX control probes
Bordet_array <- Bordet_array[!startsWith(rownames(Bordet_array), 'AFFX'),]
Bordet_RNA_tpm <- Bordet_RNA_tpm[!startsWith(rownames(Bordet_RNA_tpm), 'AFFX'),]
CIT_full <- CIT_full[!startsWith(rownames(CIT_full), 'AFFX'),]


### Quantile normalization to CIT data set + batch correction of TPM derived from kallisto
## RNA-seq
# Quantile normalization
Bordet_RNA_tpm_qnorm <- normalize.quantiles.use.target(Bordet_RNA_tpm,target = rowMeans(CIT_full))
rownames(Bordet_RNA_tpm_qnorm) <- rownames(Bordet_RNA_tpm); colnames(Bordet_RNA_tpm_qnorm) <- colnames(Bordet_RNA_tpm)

# Batch correction
batch <- ComBat(cbind(Bordet_RNA_tpm_qnorm, CIT_full),c(rep(1,ncol(Bordet_RNA_tpm_qnorm)),rep(2,ncol(CIT_full))))
Bordet_RNA_tpm_qnorm_bc <- batch[,1:ncol(Bordet_RNA_tpm_qnorm)]

## Microarray
# Quantile normalization
Bordet_array_qnorm <- normalize.quantiles.use.target(Bordet_array,target = rowMeans(CIT_full))
rownames(Bordet_array_qnorm) <- rownames(Bordet_array); colnames(Bordet_array_qnorm) <- colnames(Bordet_array)

# Batch correction
batch <- ComBat(cbind(Bordet_array_qnorm, CIT_full),c(rep(1,ncol(Bordet_array_qnorm)),rep(2,ncol(CIT_full))))
Bordet_array_qnorm_bc <- batch[,1:ncol(Bordet_array_qnorm)]


### Classification of Bordet array and RNA using the CITBCMST script
# RNA-seq
Bordet_RNA_tpm_res <- CITBCMST_modified(Bordet_RNA_tpm_qnorm_bc,Bordet_annot, plot="TRUE", title="Bordet RNA-seq samples (Qnorm + ComBat)")

# Microarray
Bordet_array_res <- CITBCMST_modified(Bordet_array_qnorm_bc,Bordet_annot, plot="TRUE", title="Bordet microarray samples (Qnorm + ComBat)")

# Comparing results
table(Bordet_array_res$classes$citbcmst.mixed==Bordet_RNA_tpm_res$classes$citbcmst.mixed)
