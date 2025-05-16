setwd("~/Dropbox (NYU Langone Health)/BGW_CJC/projects/2023_VAP/")

# Load libraries
library(dplyr)
library(lubridate)
library(tibble)
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(MatchIt)
library(edgeR)
library(tidyr)
library(decontam)
library(matrixStats)
library(forcats)
library(scales)
library(stringr)
library(writexl)
library(gtsummary)

#clear all objects
rm(list=ls(all.names=T))
gc()

# Read and load files
metadata <- read.csv(file="240326_metadata_245.csv", header=TRUE, row.names=1)
counts.airway.filter <- read.delim(file="COVIDonly_allairways_counts_taxa_postconqur.txt", header=TRUE, sep="\t", row.names = "X")
counts.host.filter <- read.delim(file="COVIDonly_RNAHost_NYUonly_01262023.txt", header=TRUE, sep="\t", row.names = "X")

# Define BASELINE (earliest lower airway sample after intubation, before VAP diagnosis)
metadata.baseline.la <- metadata %>% 
  filter(research_sample2 == "LA") %>%
  rownames_to_column('base') %>%
  group_by(study_id) %>%
  filter((vap_to_sampling < 0) %>% replace_na(TRUE)) %>% 
  filter(intub_to_sampling >= 0) %>%
  dplyr::slice(which.min(intub_to_sampling)) %>% 
  column_to_rownames('base')

# Pre-VAP (only VAP; lower airway samples after intubation, 72 hours before VAP diagnosis)
metadata.baseline.prevap <- metadata %>% 
  filter(research_sample2 == "LA") %>%
  rownames_to_column('base') %>%
  group_by(study_id) %>%
  filter(intub_to_sampling >= 2) %>%
  filter(vap_to_sampling < -3) %>%
  slice(which.max(vap_to_sampling)) %>%
  column_to_rownames('base')
metadata.baseline.prevap$prepost <- "pre_VAP"

# VAP (only VAP; lower airway samples 72 hours before & after VAP diagnosis)
metadata.baseline.vap <- metadata %>% 
  filter(research_sample2 == "LA") %>%
  rownames_to_column('base') %>%
  group_by(study_id) %>%
  filter(vap_to_sampling > -4, vap_to_sampling < 4) %>%
  slice(which.min(abs(vap_to_sampling))) %>% #changed code from which.max() to which.min(abs())
  column_to_rownames('base') # 36 samples
metadata.baseline.vap$prepost <- "VAP"

## For baseline lower airway samples between no VAP and VAP (metatranscriptome) ####
# Alpha diversity ==========================================================
nn <- intersect(rownames(metadata.baseline.la), colnames(counts.airway.filter)) # find overlapping samples
counts.airway.baseline.la <- counts.airway.filter[,nn]
counts.airway.baseline.la.rel <- sweep(counts.airway.baseline.la, 2, colSums(counts.airway.baseline.la), `/`) # relative abundance table

metadata.baseline.la$alpha_shannon <- diversity(t(counts.airway.baseline.la), index = "shannon") # calculate shannon diversity index

pairs <- list( c("VAP", "no-VAP"))
metadata.baseline.la %>% 
  mutate(had_vap = factor(had_vap, levels=c("N", "Y"), labels = c("no-VAP", "VAP"))) %>% 
  ggboxplot(x = "had_vap", y = "alpha_shannon",
          xlab = F, ylab = "Alpha diversity (Shannon)", legend = "none",
          add = "jitter", color = "had_vap", palette = c("#1874CD", "#EE6363"),
          add.params = list(size = 2, alpha = 0.5)) +
  theme(axis.line = element_line(size = 0.2)) +
  stat_compare_means(comparisons = pairs, label.y = 5, method = "wilcox.test")

# Beta diversity ==========================================================
bray_dist = vegdist(t(counts.airway.baseline.la.rel), method="bray") # calculate bray-curtis dissimilarity on relative abundance table
CmdScale <- cmdscale(bray_dist, k =10)
vars <- apply(CmdScale, 2, var)
percentVar <- round(100 * (vars/sum(vars)))
newResults <- merge(x = CmdScale, y = metadata.baseline.la, by = "row.names", all.x = TRUE)
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
centroids <- aggregate(cbind(PC1, PC2) ~ had_vap, data = newResults, mean)
newResults <- merge(newResults, centroids, by="had_vap", suffixes = c("", ".centroid"))

adonis2(bray_dist ~ metadata.baseline.la$had_vap, permutations = 999) # calculate significance

ggplot(newResults, aes(PC1, PC2, color = had_vap)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #coord_fixed() +
  scale_color_manual(values=c("#1874CD", "#EE6363")) +
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=had_vap), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=had_vap)) +
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label = c("No VAP", "VAP")), size=10) +
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background=element_blank(), axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        plot.margin=unit(c(1,1,1,1),"line"), legend.position="none") +
  annotate("text", x=0.5, y=-0.35, size=7, label="p = 0.284", fontface = 'italic')


# edgeR ==========================================================
# Differential metatranscriptome profiles between no-VAP and VAP cohorts using edgeR
# Prepare a table for the edgeR function
counts.airway.baseline.la.rename <- counts.airway.base.la
counts.airway.baseline.la.rename['had_vap',] <- metadata.baseline.la$had_vap 
counts.airway.baseline.la.rename $had_vap <- rownames(counts.airway.baseline.la.rename)
counts.airway.baseline.la.rename <- counts.airway.baseline.la.rename %>% 
  select(had_vap, everything())

counts.airway.baseline.la.newrow <- counts.airway.baseline.la.rename %>% 
  add_row(.before = 1) %>% 
  slice(-n()) %>% 
  mutate(across(everything(), ~ replace(., 1, .[n()])))
rownames(counts.airway.baseline.la.newrow)[1] <- "had_vap"
colnames(counts.airway.baseline.la.newrow) <- paste("X", seq_len(ncol(counts.airway.baseline.la.newrow)))
colnames(counts.airway.baseline.la.newrow)[1] <- "had_vap"

Edger_nonPS_KW(input_table= counts.airway.baseline.la.newrow, 
               sample_type_var_name=c("had_vap"),
               sample_types=c("N", "Y"),
               sample_type_color=c("#1874CD", "#EE6363"),  
               FDR_cut_off=0.2,
               number_display=50,
               graph_option = "volcano",
               display_all_results_volcano = "yes",
               abundance_size_volcano = "yes",
               plot_title="EdgeR: lower airway (baseline) - no VAP vs. VAP",
               output_name="EdgeR_LA_baseline_noVAP vs. VAP_volcano",
               legend_onplot="yes",
               width=15, height=10)

## For baseline lower airway samples between no VAP and VAP (host transcriptome) ####
# Alpha diversity ==========================================================
nn <- intersect(rownames(metadata.baseline.la), colnames(counts.host.filter)) # find overlapping samples
counts.host.baseline.la <- counts.host.filter[,nn]

metadata.baseline.la$alpha_shannon <- diversity(t(counts.host.baseline.la), index = "shannon") # calculate shannon diversity index

pairs <- list( c("VAP", "no-VAP"))
metadata.baseline.la %>% 
  mutate(had_vap = factor(had_vap, levels=c("N", "Y"), labels = c("no-VAP", "VAP"))) %>% 
  ggboxplot(x = "had_vap", y = "alpha_shannon",
          xlab = F, ylab = "Alpha diversity (Shannon)", legend = "none",
          add = "jitter", color = "had_vap", palette = c("#1874CD", "#EE6363"),
          add.params = list(size = 2, alpha = 0.5)) +
  theme(axis.line = element_line(size = 0.2)) +
  stat_compare_means(comparisons = pairs, label.y = 5, method = "wilcox.test")

# Beta diversity ==========================================================
bray_dist = vegdist(t(counts.host.baseline.la), method="bray")
CmdScale <- cmdscale(bray_dist, k =10)
vars <- apply(CmdScale, 2, var)
percentVar <- round(100 * (vars/sum(vars)))
newResults <- merge(x = CmdScale, y = metadata.baseline.la, by = "row.names", all.x = TRUE)
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
centroids <- aggregate(cbind(PC1, PC2) ~ had_vap, data = newResults, mean)
newResults <- merge(newResults, centroids, by="had_vap", suffixes = c("", ".centroid"))

adonis2(bray_dist ~ metadata.baseline.la$had_vap, permutations = 999) # calculate significance

ggplot(newResults, aes(PC1, PC2, color = had_vap)) +
  geom_point(size=5,alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #coord_fixed() +
  scale_color_manual(values=c("#1874CD", "#EE6363")) +
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=had_vap), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=had_vap)) +
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label = c("No VAP", "VAP")), size=10) +
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background=element_blank(), axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        plot.margin=unit(c(1,1,1,1),"line"), legend.position="none") +
  annotate("text", x=0.5, y=-0.35, size=7, label="p = 0.284", fontface = 'italic')

# edgeR ==========================================================
# Differential transcriptome profiles between no-VAP and VAP cohorts using edgeR
# Prepare a table for the edgeR function
counts.host.baseline.la.rename <- counts.host.baseline.la
counts.host.baseline.la.rename['had_vap',] <- metadata.baseline.la$had_vap 
counts.host.baseline.la.rename $had_vap <- rownames(counts.host.baseline.la.rename)
counts.host.baseline.la.rename <- counts.host.baseline.la.rename %>% 
  select(had_vap, everything())

counts.host.baseline.la.newrow <- counts.host.baseline.la.rename %>% 
  add_row(.before = 1) %>% 
  slice(-n()) %>% 
  mutate(across(everything(), ~ replace(., 1, .[n()])))
rownames(counts.host.baseline.la.newrow)[1] <- "had_vap"
colnames(counts.host.baseline.la.newrow) <- paste("X", seq_len(ncol(counts.host.baseline.la.newrow)))
colnames(counts.host.baseline.la.newrow)[1] <- "had_vap"

Edger_nonPS_KW(input_table= counts.host.baseline.la.newrow, 
               sample_type_var_name=c("had_vap"),
               sample_types=c("N", "Y"),
               sample_type_color=c("#1874CD", "#EE6363"),  
               FDR_cut_off=0.2,
               number_display=50,
               graph_option = "volcano",
               display_all_results_volcano = "yes",
               abundance_size_volcano = "yes",
               plot_title="EdgeR: lower airway (baseline) - no VAP vs. VAP",
               output_name="EdgeR_host_LA_baseline_noVAP vs. VAP_volcano",
               legend_onplot="yes",
               width=15, height=10)


