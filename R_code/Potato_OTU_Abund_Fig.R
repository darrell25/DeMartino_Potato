library(ggplot2)
library(data.table)

#This is for making the panels of the figure examining the abundance of certain OTUs in certain individuals
#This is repeated changing which CSV file is imported, the title of the figure and the name of the output file
plot_sp <- fread("T07Rfaecis.csv")

plot_order <- c("POT", "REF", "Control", "PS", "RS3")
plot_sp$Treatment <- factor(plot_sp$Treatment, plot_order)

tiff("T07Rfaecis.tiff", units="in", width=5, height=5, res=300)
ggplot(plot_sp, aes(x=Treatment, y=Rel_abun, fill = Study)) +
  geom_bar(stat = 'Identity') +
  ylab("Percent Relative Abundance") +
  ggtitle("Roseburia faecis in P21") +
  theme (plot.title = element_text(hjust = 0.5))+
  scale_fill_discrete(labels = c("in vitro", "in vivo"))
dev.off()