library(data.table)
library(dplyr)
library(lmerTest)
library(multcomp)

#This is the example for butyrate, but a simple search and replace for any other of the fermentation
#acids will generate the equivalent plot for that one (e.g. replace butyrate with propionate)


Treats_abrev <- c("BL", "Pot", "Ref")
SCFA <- c("Butyrate", "Propionate", "Acetate")

#read in metadata file that includes SCFA data
sample <- fread(file="Potato.meta.csv")
sample <- sample[order(StudyID, Diet),]

sample$Treatment <- factor(sample$Treatment, levels = Treats_abrev, labels = Treats_abrev)
sample$Treatment <- relevel(as.factor(sample$Treatment), ref = "BL")

sample.but <- sample[,c("Sample_Number", "StudyID", "Treatment", "Week", "Butyrate", "Propionate", "Acetate")]

####Treatment Effects After Averaging across weeks####

#average across weeks
sample.means <- sample.but[,lapply(.SD, mean), by=.(StudyID,Treatment), .SDcols=(c("Butyrate", "Propionate", "Acetate"))]

for (acid in SCFA) {
  #linear mixed model of butyrate values to test for differences between treatments
  input_formula <- as.formula(paste0(acid, "~ Diet + (1|StudyID)"))
  fit.scfa <- lmerTest::lmer(input_formula, data = sample.means)

  scfa_test <- glht(fit.scfa, linfct = mcp(Treatment = "Tukey"))
  
  fn <- paste0("Condition_", acid, "_results.txt")
  capture.output(summary(scfa_test, test=adjusted("BH")), file= fn, append=TRUE)
  capture.output(confint(scfa_test), file= fn, append=TRUE)
  capture.output(ls_means(fit.scfa, effect = "Diet"), file= fn, append=TRUE)
}



