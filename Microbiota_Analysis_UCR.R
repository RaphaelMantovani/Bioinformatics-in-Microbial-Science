library(ape) 
library(DESeq2) 
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(miLineage)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)

# loading and cleaning data
OTU <- read.table("C:/Users/Raphael/Downloads/jacu.final.abund.opti_mcc.0.03.norm.shared", header = TRUE, sep = '\t')
row.names(OTU) = OTU$Group
OTU.clean <- OTU[,!(names(OTU) %in% c("label", "numOtus", "Group"))]

summary <- read.table('C:/Users/Raphael/Downloads/jacu.final.abund.opti_mcc.0.03.norm.groups.summary', header = TRUE, sep ='\t')
row.names(summary) = summary$group
summary.clean <- summary[,names(summary) != c('label', 'group')]

tax <- read.table("C:/Users/Raphael/Downloads/jacu.final.abund.opti_mcc.0.03.cons.taxonomy", header = TRUE, sep = '\t')
row.names(tax) = tax$OTU
tax.clean <- tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean <- separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
tax.clean <- tax.clean[,!(names(tax.clean) %in% c('OTU', 'Size', 'Strain'))]

# ordering data so that all sets are in the same order
summary.clean <- summary.clean[order(row.names(summary.clean)),]
OTU.clean <- OTU.clean[order(row.names(OTU.clean)),]
 
# further data manipulation (if needed)
summary_test <- summary[which(summary$shannon > 0.2),]
OTU_test <- OTU.clean[which(summary$shannon > 0.2),]

# removing rare taxa
OTU.clean.abund <- OTU.clean[,which(apply(OTU.clean, 2, max) > 10)]

# calculating relative abundances
OTU.clean.relabund <- sweep(OTU.clean,1,rowSums(OTU.clean),"/")

# we can set seeds to get constant random numbers when testing data
set.seed(12891) # the number can be any that you want 

OTU.physeq = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(tax.clean))
summary.physeq = sample_data(summary.clean)
physeq.alpha = phyloseq(OTU.physeq, tax.physeq, summary.physeq)

sample_data(physeq.alpha)$shannon.physeq <- estimate_richness(physeq.alpha, measures="Shannon")
plot_richness(physeq.alpha, measures="Shannon")