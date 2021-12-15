library(phyloseq)
library(DivNet)
library(breakaway)
library(lmerTest)
library(dplyr)
library(magrittr)
library(picante)
library(data.table)
library(ggsignif)

###Data Import and Pre-processing###

OTU.matrix <- data.matrix(read.table(file="Potato_OTU.shared", header=TRUE, row.names=1, sep="\t"))
tax.matrix <- as.matrix(read.table(file="Potato_OTU.taxonomy", header=TRUE, row.names=1, sep="\t"))
sample <- read.csv(file="Potato.meta.csv", header=TRUE, row.names=1, sep=",")

#Create phyloseq object
OTU.p <- otu_table(OTU.matrix, taxa_are_rows =TRUE)
TAX.p <- tax_table(tax.matrix)
sample.p <- sample_data(sample)
tree.p <- read_tree(treefile="Potato_OTU.tree")
physeq <- phyloseq(OTU.p,TAX.p,sample.p, tree.p)

physeq_gen <- tax_glom(physeq, taxrank="Genus")

###Alpha Diverstiy Analysis via DivNet###
a_div <- physeq_gen %>% divnet(ncores = 6)
estimates <- a_div$shannon %>% summary %$% estimate
ses <- a_div$shannon %>% summary %$% error
sample$Diet <- relevel(as.factor(sample$Diet), ref="Pot")
a_div_test <- betta_random(chats = estimates,
                          ses = ses,
                          X = model.matrix(~Diet, data = sample),
                          groups=sample$StudyID)

a_div_test$table

estimates2 <- a_div$simpson %>% summary %$% estimate
ses2 <- a_div$simpson %>% summary %$% error
a_div_test2 <- betta_random(chats = estimates2,
                           ses = ses2,
                           X = model.matrix(~Diet, data = sample),
                           groups=sample$StudyID)

a_div_test2$table

###Faith's Phylogenetic Diversity Analysis via Picante###
set.seed(3488)
physeq_gen.rar <- rarefy_even_depth(physeq_gen, replace = FALSE)
OTU_gen.rar <- as.data.table(otu_table(physeq_gen.rar), keep.rownames = "OTU")
OTU_gen.rar.t <- data.table::transpose(OTU_gen.rar, keep.names = "Group", make.names = "OTU")
Faith_d.rar <- picante::pd(OTU_gen.rar.t, tree.p, include.root = FALSE)
sample$Faith <- Faith_d.rar$PD

fit.faith <- lmerTest::lmer(Faith ~ Diet + (1 | StudyID), data = sample)
anova(fit.faith)

###Plot Alpha Diversity###

#Note P-values come from tables in the results of betta_random or anova
inv_simp <- 1/estimates2

plot_sample <- data.frame(Group=rownames(sample), Diet=sample$Diet, Shannon=estimates, invSimp=inv_simp, Faith=Faith_d.rar$PD)
plot_sample$Diet <- factor(plot_sample$Diet, levels = c("BL", "Pot", "Ref") )

tiff("Genus_Shannon_Diet2.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Diet, Shannon)) + 
  geom_boxplot() + 
  ggtitle("Genus-level Shannon Diversity") + 
  theme (plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Pot", "Ref")), map_signif_level = TRUE, annotations = c("* p=0.023"), textsize = 5) +
  ylim(NA,4) +
  xlab("Condition")
dev.off()

tiff("Genus_Inverse_Simpson_Diet2.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Diet, invSimp)) + 
  geom_boxplot() +
  ggtitle("Genus-level Inverse Simpson Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Pot", "Ref")), map_signif_level = TRUE, annotations = c("** p=0.004"), textsize = 5) +
  ylim(NA,20) + 
  ylab("Inverse Simpson") +
  xlab("Condition")
dev.off()

tiff("Genus_Faith_Diet2.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sample, aes(Diet, Faith)) + 
  geom_boxplot() +
  ggtitle("Genus-level Faith's Phylogenetic Diversity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Pot", "Ref")), map_signif_level = TRUE, annotations = c("N.S. p=0.767"), textsize = 5) +
  ylim(NA,24) +
  xlab("Condition")
dev.off()

###Bray-Curtis Distance via DivNet###

#Create distance object and perform ordination
b_div_bray <- as.dist(a_div$`bray-curtis`)
b_div_bray_pcoa <- ordinate(physeq_gen, method="PCoA", distance=b_div_bray) 


tiff("DivNet_Genus_Bray_Diet2.tiff", units="in", width=5, height=5, res=300)
plot_ordination(physeq_gen, b_div_bray_pcoa, type ="samples", color = "Diet") +
  ggtitle("Bray-Curtis Dissimilarity") +
  theme (plot.title = element_text(hjust = 0.5)) +
  labs(color="Condition")
dev.off()

#Test for significance
sample.a <- as.data.frame(sample)
write("Genus Level DivNet Diet and Bray-Curtis Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = b_div_bray ~ Diet, data = sample.a), 
               file="PERMANOVA_results.txt", append=TRUE)
capture.output(adonis(formula = b_div_bray ~ Diet*Sex + Diet*Age_med + Diet*BMI_OW, 
                      data = sample.a), file="PERMANOVA_results.txt", append=TRUE)

###Aitchison Distance Analysis###

#CLR transformation of counts
physeq_gen_clr <- microbiome::transform(physeq_gen, "clr")
otu.clr <- t(as(otu_table(physeq_gen_clr), "matrix"))
#calculate aitchison distance (Euclidan distance of clr transformed data)
clr_ait <- dist(otu.clr, method='euc')
b_div_ait_pcoa <- ordinate(physeq_gen, method="PCoA", distance=clr_ait)


tiff("Genus_Aitchison_Diet2.tiff", units="in", width=5, height=5, res=300)
plot_ordination(physeq_gen, b_div_ait_pcoa, type ="samples", color = "Diet") +
  ggtitle("Aitchison Distance") +
  theme (plot.title = element_text(hjust = 0.5)) +
  labs(color="Condition")
dev.off()

write("\nGenus Level Diet and Aitchison Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = clr_ait ~ Diet, data = sample.a), 
               file="PERMANOVA_results.txt", append=TRUE)
capture.output(adonis(formula = clr_ait ~ Diet*Sex + Diet*Age_med + Diet*BMI_OW, 
                      data = sample.a), file="PERMANOVA_results.txt", append=TRUE)

###Unifrac analysis###

dist_gen_rar_uniw <- phyloseq::UniFrac(physeq_gen.rar, weighted = TRUE)
gen_rar_uniw_pcoa <- ordinate(physeq_gen.rar, method="PCoA", distance=dist_gen_rar_uniw)

tiff("Genus_UnifracW_Diet2.tiff", units="in", width=5, height=5, res=300)
plot_ordination(physeq_gen.rar, gen_rar_uniw_pcoa, type ="samples", color = "Diet") +
  ggtitle("Weighted UniFrac") +
  theme (plot.title = element_text(hjust = 0.5)) +
  labs(color="Condition")
dev.off()

write("\nGenus Level Diet and Weighted Unifrac Diversity Test", file= "PERMANOVA_results.txt", append=TRUE)
set.seed(3488)
capture.output(adonis(formula = dist_gen_rar_uniw ~ Diet, data = sample.a), 
               file="PERMANOVA_results.txt", append=TRUE)
capture.output(adonis(formula = dist_gen_rar_uniw ~ Diet*Sex + Diet*Age_med + Diet*BMI_OW, 
                      data = sample.a), file="PERMANOVA_results.txt", append=TRUE)
