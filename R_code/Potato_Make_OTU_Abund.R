library(phyloseq)
library(data.table)
library(dplyr)

OTU.df <- data.frame(read.table(file="Potato_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.df <- data.frame(read.table(file="Potato_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Potato.meta.csv", header=TRUE, row.names=1, sep=",")

##Create Species Names
tax.df$Species <- make.names (ifelse(tax.df$Genus==tax.df$Species,tax.df$Genus,paste(tax.df$Genus, tax.df$Species, sep="_")),  unique=FALSE)
tax.df$Species <- make.names (ifelse(grepl("_",tax.df$Species,fixed=TRUE ), tax.df$Species, paste(tax.df$Species, "unclassified", sep = "_")), unique = FALSE)

##Create phyloseq object##
OTU.p <- otu_table(OTU.df, taxa_are_rows =TRUE)
TAX.p <- tax_table(as.matrix(tax.df))
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

#Get the average relative abundance for each of the diet conditions 
physeq_83_5 <- transform_sample_counts(physeq_83_5, function(x) x/sum(x))

#focus on taxa of interest
taxa_list <- c("Otu00056", "Otu00098")
physeq_83_5 <- prune_taxa(taxa_list, physeq_83_5)
#focus on subject of interest
#repeat for other subjects of interest, changing the StudyID
physeq_subject <- subset_samples(physeq_83_5, StudyID == "POT 021")

physeq_BL <- data.frame (otu_table(subset_samples(physeq_subject, Diet == "BL")))
physeq_pot <- data.frame(otu_table(subset_samples(physeq_subject, Diet == "Pot")))
physeq_ref <- data.frame(otu_table(subset_samples(physeq_subject, Diet == "Ref")))

Diet.pot21 <- as.data.table(data.frame(Species=tax_table(physeq_subject)[,"Species"], BL = round(rowMeans(physeq_BL)*100,3), 
                                     POT = round(rowMeans(physeq_pot)*100,3), REF = round(rowMeans(physeq_ref)*100,3)), keep.rownames = "OTU")

fwrite(Diet.pot21, file="POT21_percent_abundances.txt", sep = "\t")