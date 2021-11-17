library(DESeq2)
library(phyloseq)

OTU.matrix <- data.matrix(read.table(file="Potato_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="Potato_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Potato.meta.csv", header=TRUE, row.names=1, sep=",")

##Create phyloseq object##
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
physeq <- phyloseq(OTU.p,TAX.p,sample.p)

###Filter sequences###

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
#gives 1,477 taxa and maintains 97.1% of reads.

###DESeq2 Analysis###

#convert to factors
sample$Diet <- factor(sample$Diet, levels = c("BL","Pot", "Ref"))

#Run analysis looking at differences in Diet
OTU83_5.matrix <- as.data.frame(otu_table(physeq_83_5))
DE_Diet_data83_5 <- DESeqDataSetFromMatrix(countData = OTU83_5.matrix, colData = sample, design = ~Diet)
DE_Diet83_5 <- DESeq(DE_Diet_data83_5)

#find the results comparing potato to refined diet
res83_5 <- results(DE_Diet83_5, contrast = c("Diet", "Pot", "Ref"), cooksCutoff = FALSE)
alpha_cut <- 0.05
sigtab83_5 <- res83_5[which(res83_5$padj < alpha_cut), ]
sigtab83_5 = cbind(as(sigtab83_5, "data.frame"), as(tax_table(physeq_83_5)[rownames(sigtab83_5), ], "matrix"))
write.csv(sigtab83_5,file="DESEQ.csv")