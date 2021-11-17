library(data.table)
library(phyloseq)
library(dplyr)


OTU.matrix <- data.matrix(read.table(file="Potato_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.dt <- fread(file = "Potato_OTU.taxonomy", header = TRUE, sep="\t")
sample <- read.csv(file="Potato.meta.csv", header=TRUE, row.names=1, sep=",")

#Make genus_species names and add OTU number

tax.dt$Species <- make.names (ifelse(tax.dt$Genus == tax.dt$Species, paste(tax.dt$Species, tax.dt$OTU, sep = "_"), 
                                paste(tax.dt$Genus, tax.dt$Species, tax.dt$OTU, sep="_")), unique=TRUE)


tax.matrix<-as.matrix(tax.dt, rownames = "OTU")

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

physeq_83 <- filter_taxa(physeq, function(x) sum(x) >=83, prune=TRUE)
prevalenceThreshold <- 0.05 * nsamples(physeq_83)

# Compute prevalence of each feature, store as data.frame
prevdf83 <- apply(X = otu_table(physeq_83),
                  MARGIN = ifelse(taxa_are_rows(physeq_83), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf83 <- data.frame(Prevalence = prevdf83,
                       TotalAbundance = taxa_sums(physeq_83),
                       tax_table(physeq_83))

keepTaxa2 <- rownames(prevdf83)[(prevdf83$Prevalence >= prevalenceThreshold)]
physeq_83_5 <- prune_taxa(keepTaxa2, physeq_83)

#Get abundances at each taxonomic rank
physeq_dom <- tax_glom(physeq_83_5, taxrank = "Domain")
physeq_phy <- tax_glom(physeq_83_5, taxrank = "Phylum")
physeq_class <- tax_glom(physeq_83_5, taxrank = "Class")
physeq_order <- tax_glom(physeq_83_5, taxrank = "Order")
physeq_fam <- tax_glom(physeq_83_5, taxrank = "Family")
physeq_gen <- tax_glom(physeq_83_5, taxrank = "Genus")

#Transform to relative abundances
physeq_dom <- transform_sample_counts(physeq_dom, function(x) x/sum(x))
physeq_phy <- transform_sample_counts(physeq_phy, function(x) x/sum(x))
physeq_class <- transform_sample_counts(physeq_class, function(x) x/sum(x))
physeq_order <- transform_sample_counts(physeq_order, function(x) x/sum(x))
physeq_fam <- transform_sample_counts(physeq_fam, function(x) x/sum(x))
physeq_gen <- transform_sample_counts(physeq_gen, function(x) x/sum(x))
physeq_sp <- transform_sample_counts(physeq_83_5, function(x) x/sum(x))

#Extract the otu and tax tables
otu_dom <- as.data.table(otu_table(physeq_dom), keep.rownames = "OTU")
tax_dom <- as.data.table(as.data.frame(tax_table(physeq_dom)), keep.rownames = "OTU")
otu_phy <- as.data.table(otu_table(physeq_phy), keep.rownames = "OTU")
tax_phy <- as.data.table(as.data.frame(tax_table(physeq_phy)), keep.rownames = "OTU")
otu_class <- as.data.table(otu_table(physeq_class), keep.rownames = "OTU")
tax_class <- as.data.table(as.data.frame(tax_table(physeq_class)), keep.rownames = "OTU")
otu_order <- as.data.table(otu_table(physeq_order), keep.rownames = "OTU")
tax_order <- as.data.table(as.data.frame(tax_table(physeq_order)), keep.rownames = "OTU")
otu_fam <- as.data.table(otu_table(physeq_fam), keep.rownames = "OTU")
tax_fam <- as.data.table(as.data.frame(tax_table(physeq_fam)), keep.rownames = "OTU")
otu_gen <- as.data.table(otu_table(physeq_gen), keep.rownames = "OTU")
tax_gen <- as.data.table(as.data.frame(tax_table(physeq_gen)), keep.rownames = "OTU")
otu_sp <- as.data.table(otu_table(physeq_sp), keep.rownames = "OTU")
tax_sp <- as.data.table(as.data.frame(tax_table(physeq_sp)), keep.rownames = "OTU")

#Convert to LEfSe formatting
Lef_file <- data.table(Sample = tax_dom$Domain)
Lef_file <- cbind(Lef_file, otu_dom[,OTU:=NULL])
phyla <- data.table(Sample = paste(tax_phy$Domain, tax_phy$Phylum, sep = "|"))
phyla <- cbind(phyla, otu_phy[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_class$Domain, tax_class$Phylum, tax_class$Class, sep = "|"))
phyla <- cbind(phyla, otu_class[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_order$Domain, tax_order$Phylum, tax_order$Class, tax_order$Order, sep = "|"))
phyla <- cbind(phyla, otu_order[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_fam$Domain, tax_fam$Phylum, tax_fam$Class, tax_fam$Order, tax_fam$Family, sep = "|"))
phyla <- cbind(phyla, otu_fam[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_gen$Domain, tax_gen$Phylum, tax_gen$Class, tax_gen$Order, tax_gen$Family, tax_gen$Genus, sep = "|"))
phyla <- cbind(phyla, otu_gen[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)
phyla <- data.table(Sample = paste(tax_sp$Domain, tax_sp$Phylum, tax_sp$Class, tax_sp$Order, tax_sp$Family, tax_sp$Genus, tax_sp$Species, sep = "|"))
phyla <- cbind(phyla, otu_sp[,OTU:=NULL])
Lef_file <- rbind(Lef_file, phyla)

sample_tb <- as.data.table(sample, keep.rownames = "Group")
sample_t <- data.table::transpose(sample_tb, make.names = "Group", keep.names = "Sample")
Lef_file <- rbind(sample_t[1:2,], Lef_file)

fwrite(Lef_file, file= "Pot_OTU_LEfSe.txt", col.name = FALSE, sep = "\t")

#version without the baseline samples
#I assume there is an easier way to do this, but this works
eval_BL <- as.logical(apply(Lef_file, 2, function(r) any(r %in% c("BL"))))
remove_BL <- Lef_file[,-(eval_BL), with = FALSE]
fwrite(remove_BL, file= "Pot_OTU_LEfSe_PR.txt", col.name = FALSE, sep = "\t")
