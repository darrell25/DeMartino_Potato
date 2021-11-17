library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(scales)

###Import data
#This is a table with each significant OTU, genus and phylum found by Ancom, DESeq2 and LEfse
#Relative abundance in each condition is included
#Significance is 1 '*' for each of the methods that detected that unit as being signficantly changed
#All significances are between the two test conditions, potato vs refined grain

da_tab <- fread("Potato_Diff_Abund.csv")
da_tab$Name <- if_else(da_tab$Level=="Species", paste(da_tab$OTU,da_tab$Name, sep="_"), da_tab$Name)

###Heatmap values
#calculate relative levels of each taxa for heatmap 
rel_nums <- da_tab[,c(4:6)]
norm_nums <- rel_nums/apply(rel_nums,1,max)

###Make the Figure
#Prepare objects for plotting
plot_vals <- cbind(Name=da_tab$Name, Level = da_tab$Level, norm_nums)
plot_sigs <- cbind(Name=da_tab$Name, Level = da_tab$Level, da_tab[,c(7:9)])
colnames(plot_sigs) <- colnames(plot_vals)

#order everything
plot_vals$Level <- factor(plot_vals$Level, levels = c("Phylum", "Genus", "Species"))
plot_vals <- plot_vals[order(plot_vals$Level, plot_vals$Name, decreasing = TRUE)]
ordernames <- plot_vals$Name

#convert to long form
plot_vals <- data.table::melt(plot_vals, variable.name = "Treatment", id.vars = c("Name", "Level"), value.name = "Rel_Change")
plot_sigs <- data.table::melt(plot_sigs, variable.name = "Treatment", id.vars = c("Name", "Level"), value.name = "Sig")

#order things so they show up in the figure correctly
plot_vals$Level <- factor(plot_vals$Level, levels = c("Phylum", "Genus", "Species"))
plot_sigs$Level <- factor(plot_sigs$Level, levels = c("Phylum", "Genus", "Species"))
treat_ord <- c("Baseline", "Potato", "Refined")
plot_vals$Treatment <- factor(plot_vals$Treatment, levels = treat_ord)
plot_sigs$Treatment <- factor(plot_sigs$Treatment, levels = treat_ord)

plot_vals <- plot_vals[order(plot_vals$Level, plot_vals$Name, plot_vals$Treatment, decreasing = TRUE)]
plot_sigs <- plot_sigs[order(plot_sigs$Level, plot_sigs$Name, plot_sigs$Treatment, decreasing = TRUE)]

plot_vals$Name <- factor(plot_vals$Name, levels = ordernames)
plot_sigs$Name <- factor(plot_sigs$Name, levels = ordernames)

#get the significant values so that they can be included on the plot
plot_sigs <- plot_sigs[is.na(plot_sigs$Sig)==FALSE,]

#make all the labels
sig_dat <- data.frame(Treatment = plot_sigs$Treatment, Name = plot_sigs$Name)
Lab_txt <- c("Baseline", "Potato", "Refined")

#make the actual plot and export it
tiff("Diff_abundance.tiff", units="in", width=7, height=6, res=300)
ggplot(plot_vals, aes(x = Treatment, y = Name)) +
  geom_tile(aes(fill = Rel_Change)) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) +
  scale_fill_gradient2(low="blue", mid = "white", high=muted("red"), midpoint = 0.5) +
  scale_x_discrete(label=Lab_txt) + 
  labs(y = "Taxa") +
  geom_text(data = sig_dat, label = plot_sigs$Sig)
dev.off()
