R Code Guide

Description of how the R code in this repository was used

Files output from the sequence processing were used along with a meta data file for the analysis in this code

Differential Abundance: Potato_OTU_ANCOM.R, Potato_OTU_DESEQ2.R and Potato_OTU_Make_LEfSe.R were used to perform differential abundance analysis by each of these methods. For LEfSe the output files were used at the Huttenhower lab Galaxy web server and run with default parameters with "StudyID" input as "Individual". Potato_Diff_Abund_Fig.R was used to then generate the figure summarizing this data. This script requires Potato_Diff_Abund.csv which summarizes and combines the results of the differential abundance methods and the relative abundances of each taxon in each condition. 

Diversity Analysis: Potato_Diversity_Genus_level.R was used to measure and compare alpha (Shannon, Inverse Simpson, Faith) and beta (Bray-Curtis, Aitchison, Weighted UniFrac) diversity across treatments and between treatments and controls.

Relative Abundance: Potato_Make_OTU_Abund.R was used to generate diet and individual specific relative abundances for certain key OTUs identified through differential abundance analysis. Potato_OTU_Abund_Fig.R was then used to generate the figures of this data. 
