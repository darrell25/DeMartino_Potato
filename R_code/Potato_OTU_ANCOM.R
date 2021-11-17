library(ggplot2)
library(nlme)
library(dplyr)
library(compositions)
library(phyloseq)
library(data.table)
library(tidyverse)
library(magrittr)

# ANCOM 2.1 from: https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R
# Place in working directory
source("ancom_v2.1.R")

OTU.matrix <- data.matrix(read.table(file="Potato_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="Potato_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Potato.meta.csv", header=TRUE, row.names=1, sep=",")

##Create phyloseq object##
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

##Create filter to have minimum 83 counts (0.001%) and 5% prevalence

physeq_83 <- filter_taxa(physeq, function(x) sum(x) >=83, prune=TRUE)

# Compute prevalence of each feature, store as data.frame
prevdf83 <- apply(X = otu_table(physeq_83),
                  MARGIN = ifelse(taxa_are_rows(physeq_83), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf83 <- data.frame(Prevalence = prevdf83,
                       TotalAbundance = taxa_sums(physeq_83),
                       tax_table(physeq_83))

prevalenceThreshold <- 0.05 * nsamples(physeq_83)
keepTaxa2 <- rownames(prevdf83)[(prevdf83$Prevalence >= prevalenceThreshold)]
physeq_83_5 <- prune_taxa(keepTaxa2, physeq_83)

#remove baseline samples
physeq_83_5 <- subset_samples(physeq_83_5, Diet != "BL")

feature_table <- as.data.frame(otu_table(physeq_83_5))
meta_data <- as.data.table(sample, keep.rownames = "SeqNum")

#setup variables for ANCOM
sample_var <- c("SeqNum")

lib_cut = 1000
neg_lb=FALSE   
main_var <- c("Diet")


# Data Pre-Processing
processed_data <- feature_table_pre_process(feature_table, meta_data, sample_var, group_var=NULL, 
out_cut = 0.05, zero_cut = 0.95, lib_cut, neg_lb)

#Run ANCOM
ANCOM_results <- ANCOM(processed_data$feature_table, processed_data$meta_data, struc_zero = processed_data$structure_zeros, main_var, p_adj_method = "BH", 
                       alpha = 0.05, adj_formula = NULL, rand_formula = "~1|StudyID")

ANCOM_results$fig

results_83_5 <- data.table(ANCOM_results$out)
results_0.6 <- results_83_5[detected_0.6==TRUE,]
result_dat <- ANCOM_results$fig$data
theresults <- filter(result_dat, taxa_id %in% results_0.6$taxa_id)
theresults <- cbind(theresults, results_0.6[,detected_0.9:detected_0.6])
theresults <- dplyr::rename(theresults, meanCLR=x)
theresults <- dplyr::rename(theresults, W_score=y)

fwrite(theresults, file = "ANCOM_OTU.sig", sep = "\t")

