# loading qiime2R

library(qiime2R)
library(RCurl)
library(phyloseq)
library(dplyr)
library(ggplot2)

table <- read_qza("C:/Users/Raphael/Downloads/table.qza")
aligned_seqs <- read_qza("C:/Users/Raphael/Downloads/masked-aligned-cl-rep-seqs.qza")
taxonomy <- parse_taxonomy(read_qza("C:/Users/Raphael/Downloads/taxonomy.qza")$data)
rooted_tree <- read_qza("C:/Users/Raphael/Downloads/rooted-tree.qza")

physeq <- qza_to_phyloseq(
  features = "C:/Users/Raphael/Downloads/table.qza",
  tree = "C:/Users/Raphael/Downloads/rooted-tree.qza", "C:/Users/Raphael/Downloads/taxonomy.qza"
)

taxonomy %>% 
  group_by(Phylum) %>% 
  filter(n() >= 3, !is.na(Phylum)) %>% 
  summarise(n_taxonomy = n()) %>% 
  arrange(desc(n_taxonomy)) %>% 
  ggplot(aes(reorder(Phylum, n_taxonomy), n_taxonomy))+
  geom_bar(stat='identity', fill = 'lightgreen', color = 'black')+
  geom_text(aes(label = n_taxonomy), size = 2.5, nudge_y = 13)+
  theme_bw()+
  ggtitle('GWAS Negative Test Phylum Count')+
  xlab(element_blank())+
  ylab(element_blank())+
  coord_flip()

taxonomy %>% 
  group_by(Family) %>% 
  filter(n() >= 5, !is.na(Family)) %>% 
  summarise(n_taxonomy = n()) %>% 
  arrange(desc(n_taxonomy)) %>% 
  head(25) %>% 
  ggplot(aes(reorder(Family, n_taxonomy), n_taxonomy))+
  geom_bar(stat='identity', fill = 'lightblue', color = 'black', width = 0.8)+
  geom_text(aes(label = n_taxonomy), size = 2.5, nudge_y = 13)+
  theme_bw()+
  ggtitle('GWAS Negative Test Family Count')+
  xlab(element_blank())+
  ylab(element_blank())+
  theme(axis.text = element_text(size = 6))+
  coord_flip()

taxonomy %>% 
  group_by(Genus) %>% 
  filter(n() >= 5, !is.na(Genus)) %>% 
  summarise(n_taxonomy = n()) %>% 
  arrange(desc(n_taxonomy)) %>% 
  head(25) %>% 
  ggplot(aes(reorder(Genus, n_taxonomy), n_taxonomy))+
  geom_bar(stat='identity', fill = 'pink', color = 'black')+
  geom_text(aes(label = n_taxonomy), size = 2.5, nudge_y = 8)+
  theme_bw()+
  ggtitle('GWAS Negative Test Genus Count')+
  xlab(element_blank())+
  ylab(element_blank())+
  theme(axis.text = element_text(size = 6))+
  coord_flip()
