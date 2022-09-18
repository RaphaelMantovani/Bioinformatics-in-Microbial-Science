setwd("C:/Users/Raphael/Documents/suen_lab/MPE-project/Final Outputs")

#Load in Libraries

library(ggplot2)

library(dplyr)

library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.

library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq

library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq

library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.

library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.

library(plotly) # A package to create interactive web graphics of use in 3D plots

library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.

library(philr) # This package provides functions for the analysis of compositional data 

library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step

library(adespatial) # Tools for the multiscale spatial analysis of multivariate data

library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks

library(qiime2R) # A package for importing qiime artifacts into an R session

library(MicrobeR) # Data visualization

library(microbiome) # Data analysis and visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(ranacapa)

library(microbiomeutilities)

library(RColorBrewer)

library(ggpubr)

library(readxl)


#Importing data

#Import the Feed Efficiency Traits into the Environment
CollapsedFeTraitsUW <- read_excel("CollapsedFeTraitsUW.xlsx")

#Create a new column containing the RFI from the most recent feed trial performed.
CollapsedFeTraitsUW$recent_RFI <- ifelse(is.na(CollapsedFeTraitsUW$rfi_3) == FALSE, CollapsedFeTraitsUW$rfi_3, ifelse(is.na(CollapsedFeTraitsUW$rfi_2) == FALSE, CollapsedFeTraitsUW$rfi_2, CollapsedFeTraitsUW$rfi_1))

#Create two values, the mean of all the RFIs and the standard deviation of the RFIs
rfi_mean <- mean(CollapsedFeTraitsUW$recent_RFI)
rfi_sd <- sd(CollapsedFeTraitsUW$recent_RFI)

#Create a new column, if the RFI value for each cow is higher than 0, then put "HE" into the new column, if it isn't (i.e. it is lower than 0) put "LE".
CollapsedFeTraitsUW$HIorLow_FE <- ifelse(CollapsedFeTraitsUW$recent_RFI > 0, "HE", "LE")

#Creating a new column, and separating the samples out into more defined groups: High Efficiency (>RFI Mean + RFI Standard Deviation), Mildly High Efficiency (>0 and <RFI Mean + RFI SD), Low Efficiency (<0-RFI Standard Deviation), and Mildly Low Efficiency (<0 and >0-RFI Standard Deviation).
CollapsedFeTraitsUW$efficiency_std <- ifelse(CollapsedFeTraitsUW$recent_RFI > (rfi_sd + rfi_mean), "HE", ifelse(CollapsedFeTraitsUW$recent_RFI >= 0 & CollapsedFeTraitsUW$recent_RFI < (rfi_sd + rfi_mean), "MHE", ifelse(CollapsedFeTraitsUW$recent_RFI < (0 - rfi_sd), "LE", "MLE")))

#Reading in the metadata.txt file
metadata = read_excel("metadata.xlsx")

#Adding new columns to the metadata table

#Adds the HE/LE categories to the metadata table. Matches the "cownumber" variable in metadata with "ID" cariable in CollapsedFeTraitsUW, and then adds the respective value from CollapsedFeTraits to metadata.
metadata$efficiency_sd = CollapsedFeTraitsUW$efficiency_std[match(metadata$cownumber,CollapsedFeTraitsUW$ID)]

#Adds the HE/MHE/LE/MLE categories to the metadata table. Follows the same method as the previous command.
metadata$HIorLow_FE = CollapsedFeTraitsUW$HIorLow_FE[match(metadata$cownumber,CollapsedFeTraitsUW$ID)]

#Adds the RFI value of the most recent trial for each animal.
metadata$recent_RFI = CollapsedFeTraitsUW$recent_RFI[match(metadata$cownumber,CollapsedFeTraitsUW$ID)]

#Adds the RFID value for each animal. RFID values are used to link to the genomic information.
metadata$RFID = CollapsedFeTraitsUW$RFID[match(metadata$cownumber,CollapsedFeTraitsUW$ID)]

#Write the metadata file to a new file in the working directory
write.table(metadata, "metadata_PlusRFI.txt", sep="\t", row.names=FALSE, quote=FALSE)


#Import the Table, rooted tree, classification file, and newly created metadata file into a phyloseq object.
physeq <- qza_to_phyloseq(features="Filtered/table-no-contam.qza",
                          tree="Phylogeny/rep_seqs_filt_aligned_masked_tree_rooted.qza", 
                          "classification.qza", 
                          metadata="metadata_PlusRFI.txt")


#Subset the Phyloseq data to have just the samples that were taken from either Arlington or Marshfield and the first sampling period.
physeq.UW_P1 <- subset_samples(physeq, (Location=="ARS" | Location=="MARS") & Sampling=="1")

#Plot alpha diversity
physeq_rarefy <- rarefy_even_depth(physeq.UW_P1, rngseed=1, sample.size=0.9*min(sample_sums(physeq.UW_P1)), replace=F)

physeq_rarefy2 <- subset_samples(physeq_rarefy, !is.na(HIorLow_FE))
  `
a.div_general <- plot_richness(physeq_rarefy2, x="HIorLow_FE", measures=c("Shannon", "simpson", "Observed"), color = "HIorLow_FE") + geom_boxplot() + xlab('Efficiency') + ggtitle('Alpha Diversity Comparison for High and Low Efficiency Samples') + scale_x_discrete(limits=c('LE','HE')) + theme_bw()
print(a.div_general)
'''
comps <- make_pairs(sample_data(physeq_rarefy2)$HIorLow_FE)

#add statistical support
a.div_general + stat_compare_means(
  comparisons = comps,
  label = "p.signif",
  tip.length = 0.05,
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
    symbols = c("xxxx", "***", "**", "*", "ns")
  ),
  method = "wilcox.test")
'''
#alpha diversity for all efficiency classifications
a.div_granular <- plot_richness(physeq_rarefy2, x="efficiency_sd", measures=c("Shannon", "simpson", "Observed"), color = "efficiency_sd") + geom_boxplot() + xlab('Efficiency') + ggtitle('Alpha Diversity Per Efficiency') + scale_x_discrete(limits=c('LE', 'MLE', 'MHE', 'HE')) + theme_bw() 
plot(a.div_granular)

#beta diversity for high and low efficiency
physeq.ord <- ordinate(physeq_rarefy2, "PCoA", "bray")
b.div.bray <- plot_ordination(physeq_rarefy2, physeq.ord, type= "samples", color= "HIorLow_FE") + geom_point(size=1.5)
b.div.bray <- b.div.bray + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic() + scale_color_brewer("Efficiency", palette = "Set2") + xlab('PC1 (24.4%)') + ylab('PC2 (6.8%)')
print(b.div.bray)

#beta diversity for all efficiency
physeq.ord <- ordinate(physeq_rarefy2, "PCoA", "bray")
b.div.bray2 <- plot_ordination(physeq_rarefy2, physeq.ord, type= "samples", color= "efficiency_sd") + geom_point(size=1.5)
b.div.bray2 <- b.div.bray2 + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic() + scale_color_brewer("Efficiency", palette = "Set2") + xlab('PC1 (24.4%)') + ylab('PC2 (6.8%)')
print(b.div.bray2)
'''
# convert to relative abundance
physeq_rel <- microbiome::transform(physeq, "compositional")
physeq.ord.wuni <- ordinate(physeq_rel, "PCoA", "unifrac", weighted=T)
b.div.wuni <- plot_ordination(physeq_rel, physeq.ord.wuni, type= "samples", color= "HIorLow_FE") + geom_point(size=1.5)
b.div.wuni <- b.div.wuni + stat_ellipse() + ggtitle("Weighted Unifrac")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.wuni)
'''
# convert phyloseq to deseq
bsdds <- phyloseq_to_deseq2(physeq_rarefy2, ~ HIorLow_FE)
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(bsdds), 1, gm_mean)
bsdds <- estimateSizeFactors(bsdds, geoMeans = geoMeans)
# DeSeq function tests for differential abundance 
bsdds <- DESeq(bsdds, test="Wald", fitType="parametric")
# Results function call creates a table of the results of the tests
res <- results(bsdds, cooksCutoff = FALSE)
sigtab <- res[which(res$padj < 1), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_rarefy2)[rownames(sigtab), ], "matrix"))

# Cleaning up the table a little for legibility
posigtab <- sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab <- posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
# Bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands
sigtabgen <- subset(sigtab, !is.na(Genus))
# Phylum order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))