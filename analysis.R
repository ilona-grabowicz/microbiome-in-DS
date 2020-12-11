### Figures:
# Differentially abundant species (DA) - DeSeq2
# Plotting the top DA species 
# Plotting the most abundant species in DS and in WT mice
# Venn diagrams of species present in WT and DS
# Multidimensional Scaling Analysis (MDS) 
# Hierarchical Clustering (HCL)
# Correlation of species within WT and DS mice
# Validation with total RNA HiSeq data - presence of species - Venn diagram
# Correlation of species abundances in 16S and total RNA HiSeq data

### Supplementary Data and Figures:
# Numbers of reads per sample in 16S data
# Numbers of reads per sample in total RNA data 
# Rarefaction curves in 16S data 
# Rarefaction curves in total RNA data 
# Validation with total RNA HiSeq data - presence of species - Venn diagram for each condition
# Correlation of species within WT and DS mice - validation in total RNA data

### Comparison with metadata and/or with the human data:
# Compare numbers of detected species in mice and human data (Biagi et al. 2014) (Venn diagram)
# Compare numbers of species (mice and humans) per genotype
# Correlate mouse with mouse abundances data and human with human, and compare them
# Clustering of the human data alone

# ANALYSIS:
#################### Differential taxa abundance in 16S data analysis - with DeSeq2 ####################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
library("DESeq2")
# Species level
cts <- read.csv('split_file_combined_16S_S.csv', header=T, stringsAsFactors = F, row.names = 1) # choose a taxonomic level: 'P', 'C', 'O', 'F', 'G', 'S'
cts$level<-NULL # Counts table
cts$X<- NULL
colnames(cts) <- sub("X", "", colnames(cts))
coldata <- data.frame(matrix(nrow=ncol(cts),ncol=2)) # Metadata table
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c('strain', 'day')
coldata$day <- factor(rep(c('4', '0'),12))
coldata$strain <- factor(c(rep('WT', 12), rep('DS', 12)))

dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata,
                              design = ~ strain)
keep <- rowMeans(counts(dds)) > 20 # removing taxa with low counts
dds <- dds[keep,]
#dds$day <- factor(dds$day, levels = c("0","4")) # timepoints analysis
dds$strain <- factor(dds$strain, levels = c("WT","DS")) # mouse strains analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered_S <- data.frame(res[order(res$pvalue),])
write.csv(resOrdered_S, 'S_genotype.csv')

# Phylum level:
cts <- read.csv('split_file_combined_16S_P.csv', header=T, stringsAsFactors = F, row.names = 1) # choose a taxonomic level: 'P', 'C', 'O', 'F', 'G', 'S'
cts$level<-NULL # Counts table
cts$X<- NULL
colnames(cts) <- sub("X", "", colnames(cts))
coldata <- data.frame(matrix(nrow=ncol(cts),ncol=2)) # Metadata table
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c('strain', 'day')
coldata$day <- factor(rep(c('4', '0'),12))
coldata$strain <- factor(c(rep('WT', 12), rep('DS', 12)))

dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata,
                              design = ~ strain)
keep <- rowMeans(counts(dds)) > 20 # removing taxa with low counts
dds <- dds[keep,]
#dds$day <- factor(dds$day, levels = c("0","4")) # timepoints analysis
dds$strain <- factor(dds$strain, levels = c("WT","DS")) # mouse strains analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered_P <- data.frame(res[order(res$pvalue),])
write.csv(resOrdered_P, 'P_genotype.csv')

### Comparing diets separately in WT and DS mice
# WT:
cts <- read.csv('split_file_combined_16S_S.csv', header=T, stringsAsFactors = F, row.names = 1) # choose a taxonomic level: 'P', 'C', 'O', 'F', 'G', 'S'
cts$level<-NULL # Counts table
cts$X<- NULL
colnames(cts) <- sub("X", "", colnames(cts))
cts_WT <- cts[,grepl('6684|6686|6687|6690|6691|6694', colnames(cts))]
coldata_WT <- data.frame(matrix(nrow=ncol(cts_WT),ncol=1)) # Metadata table
rownames(coldata_WT) <- colnames(cts_WT)
colnames(coldata_WT) <- c('day')
coldata_WT$day <- factor(rep(c('4', '0'),6))

dds <- DESeqDataSetFromMatrix(countData = cts_WT, 
                              colData = coldata_WT,
                              design = ~ day)
keep <- rowMeans(counts(dds)) > 20 # removing taxa with low counts
dds <- dds[keep,]
dds$day <- factor(dds$day, levels = c("0","4")) # timepoints analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered_WT <- data.frame(res[order(res$pvalue),])
write.csv(resOrdered_WT, 'S_diet_WT.csv')

# DS:
cts_DS <- cts[,grepl('6696|6697|6698|6700|6701|6702', colnames(cts))]
coldata_DS <- data.frame(matrix(nrow=ncol(cts_DS),ncol=1)) # Metadata table
rownames(coldata_DS) <- colnames(cts_DS)
colnames(coldata_DS) <- c('day')
coldata_DS$day <- factor(rep(c('4', '0'),6))

dds <- DESeqDataSetFromMatrix(countData = cts_DS, 
                              colData = coldata_DS,
                              design = ~ day)
keep <- rowMeans(counts(dds)) > 20 # removing taxa with low counts
dds <- dds[keep,]
dds$day <- factor(dds$day, levels = c("0","4")) # timepoints analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered_DS <- data.frame(res[order(res$pvalue),])
write.csv(resOrdered_DS, 'S_diet_DS.csv')


#################### Plotting the top DA species #################### 
library(plyr)
library(data.table)
library(ggplot2)
library(scales) 
library(dplyr)
data <- read.csv('split_file_combined_normalized_16S_S.csv', header=T, stringsAsFactors = F) 
data$level <- data$X <- NULL
colnames(data) <- c('taxon', 'WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')

DA_species <- row.names(resOrdered_S[2,])
#DA_species <- row.names(resOrdered[21,])
data_plot_DA <- data[data$taxon==DA_species,]
data_plot_DA <- melt(data_plot_DA)
data_plot_DA$strain <- NA
data_plot_DA$strain[ grepl('WT', data_plot_DA$variable)] <- 'WT'
data_plot_DA$strain[ grepl('DS', data_plot_DA$variable)] <- 'Ts65Dn'

ggplot(data_plot_DA, aes(x = strain, y = value, color = strain)) +
  geom_jitter(size=5, alpha = 0.5, width=0.2) +
  scale_y_continuous(labels=function(x) format(x, big.mark = " ", scientific = FALSE)) +
  guides(color=guide_legend(title="Mouse genotype")) +
  xlab('Mouse genotype') +
  ylab(paste0(DA_species, '\nnormalized counts'))
ggsave(paste0(DA_species, '_normalized_counts.png'))
ggsave(paste0(DA_species, '_normalized_counts.svg'))

  #################### Plotting the most abundant species in DS and in WT mice #################### 

## Species ##
data <- read.csv('split_file_combined_normalized_16S_S.csv', header=T, stringsAsFactors = F) 
data$level <- data$X <- NULL
colnames(data) <- c('taxon', 'WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')

data_plot<- data
data_plot$median <- apply(data_plot[grepl('_', colnames(data_plot))], 1, median)
data_plot <- arrange(data_plot, desc(median))
data_plot <- data_plot[1:20,]
data_plot$median <- NULL
taxon_order <- c(data_plot$taxon)
data_plot <- melt(data_plot)
data_plot$strain <- NA
data_plot$strain[ grepl('WT', data_plot$variable)] <- 'WT'
data_plot$strain[ grepl('DS', data_plot$variable)] <- 'Ts65Dn'

data_plot$taxon <- factor(as.character(data_plot$taxon), levels = taxon_order) # setting order of the x axis
data_plot$significance <- NA
for (i in 1:nrow(data_plot)) {
  species <- as.character(data_plot$taxon[i])
  print(species)
  if (!is.na(resOrdered_S[species, "padj"])) {
    if (resOrdered_S[species, "padj"]<0.05) {
    data_plot$significance[i] <- '*'
    }
  }
}

ggplot(data_plot, aes(x=taxon, y=value, fill=strain)) + 
   geom_boxplot() +
   scale_y_continuous(trans='log10',
                      breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
   xlab('Species') +
   ylab( 'Log10(normalized mean number of \nreads / sample)') +
   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12)) +
  guides(fill=guide_legend(title="Mouse strain")) +
  geom_text(x = 1:nrow(data_plot), y = 5.2, label=data_plot$significance, size=7) 
ggsave('Top_20_abundant_species.svg')

## Phyla ##
data_phyla <- read.csv('split_file_combined_normalized_16S_P.csv', header=T, stringsAsFactors = F) 
data_phyla$level <- data_phyla$X <- NULL
colnames(data_phyla) <- c('taxon', 'WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')

data_plot<- data_phyla
data_plot$median <- apply(data_plot[grepl('_', colnames(data_plot))], 1, median)
data_plot <- arrange(data_plot, desc(median))
data_plot <- data_plot[1:20,]
data_plot$median <- NULL
taxon_order <- c(data_plot$taxon)
data_plot <- melt(data_plot)
data_plot$strain <- NA
data_plot$strain[ grepl('WT', data_plot$variable)] <- 'WT'
data_plot$strain[ grepl('DS', data_plot$variable)] <- 'Ts65Dn'

data_plot$taxon <- factor(as.character(data_plot$taxon), levels = taxon_order) # setting order of the x axis
data_plot$significance <- NA
for (i in 1:nrow(data_plot)) {
  species <- as.character(data_plot$taxon[i])
  print(species)
  if (!is.na(resOrdered_P[species, "padj"])) {
    if (resOrdered_P[species, "padj"]<0.05) {
      data_plot$significance[i] <- '*'
    }
  }
}

ggplot(data_plot, aes(x=taxon, y=value, fill=strain)) + 
  geom_boxplot() +
  scale_y_continuous(trans='log10',
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab('Phyla') +
  ylab( 'Log10(mean number of reads / sample)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12)) +
  guides(fill=guide_legend(title="Mouse strain")) +
  geom_text(x = 1:nrow(data_plot), y = 5.2, label=data_plot$significance, size=7)
ggsave('Top_20_abundant_phyla.svg')

#################### Venn diagrams #################### 
library(venn)

data_venn <- data[,c(2:ncol(data))]
rownames(data_venn) <- data$taxon
data_venn <- data_venn[rowSums(data_venn)>19,] #removing low-count species
DS <- row.names(data_venn[rowSums(data_venn[,grepl('DS.*T0', colnames(data_venn))])>0,])
WT <- row.names(data_venn[rowSums(data_venn[,grepl('WT.*T0', colnames(data_venn))])>0,])

input <- list(DS, WT)
pdf('venn_WT_DS_control.pdf', width = 4, height = 4)
venn(input, snames = c('DS', 'WT'), zcolor = c('#F8766D', '#00BFC4'))
dev.off()

#library(eulerr)
#plot(euler(list(DS = DS, WT = WT)), fills = c('#F8766D', '#00BFC4'))

#################### Multidimensional Scaling Analysis #################### 

library(vegan)
data_mds <- data
rownames(data_mds) <- data_mds$taxon
data_mds$taxon <- NULL
data_mds <- data_mds[rowSums(data_mds)>19,] #removing low-count species
data_mds<-data.frame(t(data_mds))

data.classes.genotype<-c(rep('#00BFC4', 12), rep('#F8766D', 12)) # c(rep('DN', 12), rep('WT', 12))
data.classes.day<-c(rep(c(1, 2),12)) # c(rep(c('T0', 'T4'),12))

# https://websites.pmc.ucsc.edu/~mclapham/Rtips/ordination.htm
PCoA.res<-capscale(data_mds~1,distance="canberra", dfun=vegdist) #bray
eig<-PCoA.res$CA$eig # percentages of each MDS component
percentages_explained <- eig / sum(eig)

png('MDS_species_Bray_strain_MDS3.png')
plot(scores(PCoA.res, choices=c(3,4), display="sites"), col=data.classes.genotype, pch=19, xlab = paste(c("MDS3 \n", as.character(round(percentages_explained[[3]]*100, 2)) , "%"), collapse=""))
#dis <- vegdist(data_mds,method="canberra")
#text(scores(PCoA.res, choices=c(3,4), display="sites"), labels(dis), col=data.classes.genotype, pos=1)
mtext(side = 2, text = paste(c(as.character(round(percentages_explained[[4]]*100, 2)) , "%"), collapse=""), line = 2)
legend("topleft", c('WT', 'Ts65Dn'), pch=19, col=c('#00BFC4','#F8766D'))
dev.off()

scores<-data.frame(scores(PCoA.res, choices=c(3,4), display="species"))
order<-order(abs(scores$MDS3), decreasing=TRUE)
scores[order,]
write.table(scores[order,], file="Scores_MDS_Bray_S_strains.csv", sep=',', row.names = TRUE, col.names = TRUE)


#################### Hierarchical Clustering #################### 

library(vegan)
library(dendextend)

data_hcl <- data
rownames(data_hcl) <- data_hcl$taxon
data_hcl$taxon <- NULL
data_hcl <- data_hcl[rowSums(data_hcl)>19,]
data_hcl<-data.frame(t(data_hcl))
dist_matr <- vegdist(data_hcl,method="bray") 
dend<-as.dendrogram(hclust(dist_matr, method='ward.D')) ##
labels_colors(dend)<-c(rep('#00BFC4', 8), '#F8766D', rep('#00BFC4', 3), '#F8766D', '#00BFC4', rep('#F8766D', 10))
png(file="HCL_strain_S.png")
par(mar=c(5,4.1,2.1,2.1))
plot(dend)
legend("topright", c('WT', 'Ts65Dn'), pch=19, col=c('#00BFC4', '#F8766D'))
dev.off()

dev.off() # Resetting par()

#################### Correlation of species within WT and DS mice #################### 

library(ggplot2)
library(scales)
library(ggrepel)
library(data.table)
library(reshape2)
library(grid)

data_corr <- data
row.names(data_corr) <- data_corr$taxon
data_corr$taxon <- NULL
data_corr <- data_corr[rowSums(data_corr)>19,]

WT_T0 <- colnames(data_corr)[grepl('WT.*T0', colnames(data_corr))]
mice1 <- WT_T0
mice2 <- mice1 [ mice1 != mice1[1]]
WT_T0_kor <- c()
WT_T0_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr[,i], data_corr[,j])
      WT_T0_kor <- c(WT_T0_kor, k$estimate)
      WT_T0_p <- c(WT_T0_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

WT_T4 <- colnames(data_corr)[grepl('WT.*T4', colnames(data_corr))]
mice1<- WT_T4
mice2 <- mice1 [ mice1 != mice1[1]]
WT_T4_kor <- c()
WT_T4_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr[,i], data_corr[,j])
      WT_T4_kor <- c(WT_T4_kor, k$estimate)
      WT_T4_p <- c(WT_T4_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

DS_T0 <- colnames(data_corr)[grepl('DS.*T0', colnames(data_corr))]
mice1<- DS_T0
mice2 <- mice1 [ mice1 != mice1[1]]
DS_T0_kor <- c()
DS_T0_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr[,i], data_corr[,j])
      DS_T0_kor <- c(DS_T0_kor, k$estimate)
      DS_T0_p <- c(DS_T0_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

DS_T4 <- colnames(data_corr)[grepl('DS.*T4', colnames(data_corr))]
mice1<- DS_T4
mice2 <- mice1 [ mice1 != mice1[1]]
DS_T4_kor <- c()
DS_T4_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr[,i], data_corr[,j])
      DS_T4_kor <- c(DS_T4_kor, k$estimate)
      DS_T4_p <- c(DS_T4_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}
# 
# WT_T0_kor
# WT_T4_kor
# DS_T0_kor
# DS_T4_kor

data_plot_corr <- data.frame('WT_T0' = c(WT_T0_kor),
                   'WT_T4' = c(WT_T4_kor),
                   'Ts65Dn_T0' = c(DS_T0_kor),
                   'Ts65Dn_T4' = c(DS_T4_kor) )

data_plot_corr<-melt(data_plot_corr)

ggplot(data_plot_corr, aes(x=variable, y=value)) +
  geom_boxplot(fill = c('#00BFC4', '#00BFC4', '#F8766D', '#F8766D'), alpha =0.8) +
  ylab('Pearson R') + xlab('Condition') +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_x_discrete (labels = c('WT T0', 'WT T4', 'Ts65Dn T0', 'Ts65Dn T4'))# +
 # geom_signif(comparisons = list(c("Ts65Dn_T0", "Ts65Dn_T4"),
 #                                c('WT_T0', 'WT_T4')), 
#              map_signif_level=TRUE)

ggsave('correlations_of_samples_within_conditions1.png')

wilcox.test(c(WT_T0_kor, WT_T4_kor), c(DS_T0_kor, DS_T4_kor))
wilcox.test(c(WT_T0_kor, DS_T0_kor), c(WT_T4_kor, DS_T4_kor))
median(c(as.numeric(WT_T0_kor), as.numeric(WT_T4_kor)))
median(c(as.numeric(WT_T0_p), as.numeric(WT_T4_p)))
median(c(as.numeric(DS_T0_kor), as.numeric(DS_T4_kor)))
median(c(as.numeric(DS_T0_p), as.numeric(DS_T4_p)))
median(c(as.numeric(DS_T0_kor)))
median(c(as.numeric(DS_T4_kor)))
median(c(as.numeric(WT_T0_kor)))
median(c(as.numeric(WT_T4_kor)))

################### Validation with total RNA HiSeq data - presence of species ################### 

library(venn)

data_venn <- data[,c(2:ncol(data))]
rownames(data_venn) <- data$taxon
#data_venn <- data_venn[rowSums(data_venn)>19,] #removing low-count species

data_total_rna <- read.csv('split_file_combined_total_RNA_S.csv', header=T, stringsAsFactors = F, row.names = 1) 
data_total_rna$level <- NULL
colnames(data_total_rna) <- c('WT7_T3', 'DS7_T3', 'DS7_T1', 'DS8_T3', 'WT8_T3', 'WT4_T3', 'WT4_T1', 'WT7_T1', 'DS8_T1')

sequencing_16S <- row.names(data_venn)
sequencing_total_RNA <- row.names(data_total_rna)
input <- list(sequencing_16S, sequencing_total_RNA)
# pdf('venn_sequencing_16S_and_totalRNA.pdf', width = 4, height = 4)
# detach("package:eulerr", unload=TRUE)
# venn(input, snames = c('16S', 'total RNA'), zcolor = c('green', 'gray'))
# dev.off()

library(eulerr)
pdf('venn_sequencing_16S_and_totalRNA.pdf', width = 4, height = 4)
plot(euler(list('16S' = sequencing_16S, 'total RNA' = sequencing_total_RNA)), fills = c('green', 'gray'))
dev.off()

################### Correlation of species abundances in 16S and total RNA HiSeq data ################### 

library(ggplot2)
library(scales)
library(ggrepel)
library(grid)

data_16s <- data
row.names(data_16s) <- data_16s$taxon
data_16s$taxon <- NULL
data_16s <- data_16s[rowSums(data_16s)>19,]
data_16s$mean <- apply(data_16s, 1, mean)

data_total_rna <- read.csv('split_file_combined_total_RNA_S.csv', header=T, stringsAsFactors = F, row.names = 1) 
colnames(data_total_rna) <- c('level', 'WT7_T3', 'DS7_T3', 'DS7_T1', 'DS8_T3', 'WT8_T3', 'WT4_T3', 'WT4_T1', 'WT7_T1', 'DS8_T1', 'X')
data_total_rna$level <- data_total_rna$X <- NULL
data_total_rna <- data_total_rna[rowSums(data_total_rna)>10,]
data_total_rna$mean <- apply(data_total_rna, 1, mean)

species <- rownames(data_total_rna)[rownames(data_total_rna)%in%rownames(data_16s)]

b2 <- b3 <-c()
for (sp in species) {
  b2 <- c(b2, data_total_rna$mean[rownames(data_total_rna)==sp])
  b3 <- c(b3, data_16s$mean[rownames(data_16s)==sp])
}

korelacja <- cor.test(b2, b3)
grob <- grobTree(textGrob(paste0('Pearson R = ', round(korelacja$estimate, digits=2), '\np-val = ', signif(korelacja$p.value, digits=3) ), x=0.1,  y=0.95, hjust=0, gp=gpar(fontface="bold")))
data <- data.frame(b2, b3)
ggplot(data, aes(x=b2, y=b3)) +
  geom_point() +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)), name='Total RNA: log2 (mean of species abundances in all mice)') +
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)), name='16S: log2 (mean of species abundances in all mice)') +
  annotation_custom(grob) 
ggsave('correlation_of_data_total_rna_and_data_16s_species_abundances_means_over_mice.png')



####### SUPPLEMENTARY FIGURES / NUMBERS ####### :

#################### Numbers of reads per sample in 16S data #################### 

all_counts <- read.csv('files_combined_SK_16S.csv', stringsAsFactors = T, header=T, nrows = 2, row.names = 1)
all_counts$X <- all_counts$level <- NULL
all_reads <- all_counts["root",] +all_counts["unclassified",] 
sum(all_reads)
median(as.numeric(all_reads))
min(as.numeric(all_reads))
max(as.numeric(all_reads))

#################### Numbers of reads per sample in total RNA data #################### 

all_counts_total <- read.csv('files_combined_SK_total_paired.csv', header=T, nrow = 2, stringsAsFactors = F, row.names = 1) 
all_counts_total$X <- all_counts_total$level <- NULL
all_reads_total <- all_counts_total["root",] +all_counts_total["unclassified",] 
sum(all_reads_total)
median(as.numeric(all_reads_total))
min(as.numeric(all_reads_total))
max(as.numeric(all_reads_total))

################### Rarefaction curves in 16S data #################### 

library(randomcoloR)
data_rarefaction <- read.csv('split_file_combined_16S_S.csv', header=T, stringsAsFactors = F, row.names = 1)
data_rarefaction$level <- data_rarefaction$X <- NULL
colnames(data_rarefaction) <- c('WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')


funkcja <- function (repeats, sample_min, sample_max, by_how_many) {
  counting <- means <- c()
  x <- seq(sample_min, sample_max, by=by_how_many)
  for (sampl in x) {
    counting <- c()
    for (rep in 1:repeats) {
      species_nr <- length(unique(sample(species, size = sampl, prob = probs, replace=T)))
      counting <- c(counting, species_nr)
    }
    means <- c(means, mean(counting))
  }
  y <- means
  return(list(sample_max,y))
}
species <- rownames(data_rarefaction)

svg('Rarefaction_curve_with_species_16S.svg')
plot (seq(0, 250000, by=10000), c(rep(0,25),400), type = 'n', xlab = 'Number of sequences drawn', ylab='Number of unique species', main='Rarefaction curve for species')
for (kolumna in 1:ncol(data_rarefaction)) {
  print(colnames(data_rarefaction[kolumna]) )
  probs <- data_rarefaction[,kolumna] / sum(data_rarefaction[,kolumna])
  serie <- funkcja (20, 0, sum(data_rarefaction[,kolumna]), 10000)[[2]]
  lines(seq(0, sum(data_rarefaction[,kolumna]), by=10000), serie, col=randomColor(kolumna))
}
#legend('right', leg.txt, col=colors, lty=c(1,1))
dev.off()

################### Rarefaction curves in total RNA data #################### 

library(randomcoloR)
data_rarefaction <- read.csv('split_file_combined_total_RNA_S.csv', header=T, stringsAsFactors = F, row.names = 1)
colnames(data_rarefaction) <- c('level', 'WT7_T3', 'DS7_T3', 'DS7_T1', 'DS8_T3', 'WT8_T3', 'WT4_T3', 'WT4_T1', 'WT7_T1', 'DS8_T1', 'X')
data_rarefaction$level <- data_rarefaction$X <- NULL

funkcja <- function (repeats, sample_min, sample_max, by_how_many) {
  counting <- means <- c()
  x <- seq(sample_min, sample_max, by=by_how_many)
  for (sampl in x) {
    counting <- c()
    for (rep in 1:repeats) {
      species_nr <- length(unique(sample(species, size = sampl, prob = probs, replace=T)))
      counting <- c(counting, species_nr)
    }
    means <- c(means, mean(counting))
  }
  y <- means
  return(list(sample_max,y))
}
species <- rownames(data_rarefaction)

svg('Rarefaction_curve_with_species_total_RNA.svg')
plot (seq(0, 10000000, by=1000000), c(rep(0,10),2000), type = 'n', xlab = 'Number of sequences drawn', ylab='Number of unique species', main='Rarefaction curve for species')
for (kolumna in 1:ncol(data_rarefaction)) {
  print(colnames(data_rarefaction[kolumna]) )
  probs <- data_rarefaction[,kolumna] / sum(data_rarefaction[,kolumna])
  serie <- funkcja (1, 0, sum(data_rarefaction[,kolumna]), 500000)[[2]]
  lines(seq(0, sum(data_rarefaction[,kolumna]), by=500000), serie, col=randomColor(kolumna))
}
#legend('right', leg.txt, col=colors, lty=c(1,1))
dev.off()


#################### Validation with total RNA HiSeq data - presence of species - Venn diagram for each condition #################### 

library(venn)
data <- read.csv('split_file_combined_normalized_16S_S.csv', header=T, stringsAsFactors = F) 
data$level <- data$X <- NULL
colnames(data) <- c('taxon', 'WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')

data_venn <- data[,c(2:ncol(data))]
rownames(data_venn) <- data$taxon

data_total_rna <- read.csv('split_file_combined_total_RNA_S.csv', header=T, stringsAsFactors = F, row.names = 1) 
colnames(data_total_rna) <- c('level', 'WT7_T3', 'DS7_T3', 'DS7_T1', 'DS8_T3', 'WT8_T3', 'WT4_T3', 'WT4_T1', 'WT7_T1', 'DS8_T1', 'X')
data_total_rna$level <- data_total_rna$X <- NULL

# WT
sequencing_16S <- data_venn[,colnames(data_venn[grepl('WT', colnames(data_venn))])]
sequencing_16S <- sequencing_16S[rowSums(sequencing_16S)>0,]
sequencing_16S <- rownames(sequencing_16S)

sequencing_total_RNA <- data_total_rna[,colnames(data_total_rna[grepl('WT', colnames(data_total_rna))])]
sequencing_total_RNA <- sequencing_total_RNA[rowSums(sequencing_total_RNA)>0,]
sequencing_total_RNA <- row.names(sequencing_total_RNA)

input <- list(sequencing_16S, sequencing_total_RNA)


library(eulerr)
pdf('venn_sequencing_16S_and_totalRNA_WT.pdf', width = 4, height = 4)
plot(euler(list('16S' = sequencing_16S, 'total RNA' = sequencing_total_RNA)), fills = c('green', 'gray'))
dev.off()

# detach("package:eulerr", unload=TRUE)
# pdf('venn_sequencing_16S_and_totalRNA_WT.pdf', width = 4, height = 4)
# venn(input, snames = c('16S', 'total RNA'), zcolor = c('green', 'gray'))
# dev.off()

# DS

sequencing_16S <- data_venn[,colnames(data_venn[grepl('DS', colnames(data_venn))])]
sequencing_16S <- sequencing_16S[rowSums(sequencing_16S)>0,]
sequencing_16S <- rownames(sequencing_16S)

sequencing_total_RNA <- data_total_rna[,colnames(data_total_rna[grepl('DS', colnames(data_total_rna))])]
sequencing_total_RNA <- sequencing_total_RNA[rowSums(sequencing_total_RNA)>0,]
sequencing_total_RNA <- row.names(sequencing_total_RNA)

input <- list(sequencing_16S, sequencing_total_RNA)

library(eulerr)
pdf('venn_sequencing_16S_and_totalRNA_DS.pdf', width = 4, height = 4)
plot(euler(list('16S' = sequencing_16S, 'total RNA' = sequencing_total_RNA)), fills = c('green', 'gray'))
dev.off()
# 
# detach("package:eulerr", unload=TRUE)
# pdf('venn_sequencing_16S_and_totalRNA_DS.pdf', width = 4, height = 4)
# venn(input, snames = c('16S', 'total RNA'), zcolor = c('green', 'gray'))
# dev.off()

#################### Correlation of species within WT and DS mice - validation in total RNA HiSeq data #################### 

library(ggplot2)
library(scales)
library(ggrepel)
library(data.table)
library(reshape2)
library(grid)

data_total_rna <- data_total_rna
data_total_rna <- data_total_rna[rowSums(data_total_rna)>19,]

# Only WT and DS

WT_total <- colnames(data_total_rna)[grepl('WT', colnames(data_total_rna))]
mice1 <- WT_total
mice2 <- mice1 [ mice1 != mice1[1]]
WT_total_kor <- c()
WT_total_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_total_rna[,i], data_total_rna[,j])
      WT_total_kor <- c(WT_total_kor, k$estimate)
      WT_total_p <- c(WT_total_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

DS_total <- colnames(data_total_rna)[grepl('DS', colnames(data_total_rna))]
mice1<- DS_total
mice2 <- mice1 [ mice1 != mice1[1]]
DS_total_kor <- c()
DS_total_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_total_rna[,i], data_total_rna[,j])
      DS_total_kor <- c(DS_total_kor, k$estimate)
      DS_total_p <- c(DS_total_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

data_plot_corr_total <- setNames(data.frame(matrix(ncol = 2, nrow = length(WT_total_kor))), c("WT_total", "DS_total"))
data_plot_corr_total$WT_total <- WT_total_kor
data_plot_corr_total$DS_total[1:length(DS_total_kor)] <- DS_total_kor

data_plot_corr_total<-melt(data_plot_corr_total)

ggplot(data_plot_corr_total, aes(x=variable, y=value)) +
  geom_boxplot(fill = c('#00BFC4', '#F8766D'), alpha =0.8) +
  ylab('Pearson R') + xlab('Mouse strain') +
  scale_x_discrete (labels = c('WT', 'DS')) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))

ggsave('correlations_of_samples_within_conditions_total_RNA.png')

wilcox.test(WT_total_kor, DS_total_kor)

median(as.numeric(WT_total_kor))
median(as.numeric(DS_total_kor))

median(as.numeric(WT_total_p))
median(as.numeric(DS_total_p))

#################### Comparison with metadata and/or with the human data: ####################


#################### Compare numbers of detected species in mice and human data (Biagi et al. 2014) (Venn diagram) #################### 

# Numbers of species (Venn)
human <- read.csv('./split_file_combined_Biagi_S.csv', stringsAsFactors = F, header=T)
human$level<-human$X <- NULL
colnames(human)<-c('taxon', paste0(rep('DOWN0', 9), 1:9), paste0(rep('DOWN', 8), 10:17), paste0(rep('CTRL0', 9), 1:9),  paste0(rep('CTRL', 7), 10:16))
human <- human[rowSums(human[,2:ncol(human)])>19,]

data <- read.csv('split_file_combined_normalized_16S_S.csv', header=T, stringsAsFactors = F) 
data$level <- data$X <- NULL
data <- data[rowSums(data[,2:ncol(data)])>19,]

library(eulerr)
svg('venn_human_mice.svg')
plot(euler(list('Human\nn=142' = human$taxon, 'Mice\nn=353' = data$taxon)), fills = c('#7CAE00', '#C77CFF'))
dev.off()

#################### Compare numbers of species (mice and humans) per genotype ####################
library(eulerr)
human <- read.csv('./split_file_combined_normalized_Biagi_S.csv', stringsAsFactors = F, header=T)
human$level<-NULL
colnames(human)<-c('taxon', paste0(rep('DOWN0', 9), 1:9), paste0(rep('DOWN', 8), 10:17), paste0(rep('CTRL0', 9), 1:9),  paste0(rep('CTRL', 7), 10:16))

data <- read.csv('split_file_combined_normalized_16S_S.csv', header=T, stringsAsFactors = F) 
data$level <- data$X <- NULL

human_DS <- human[,grepl('DOWN|taxon', colnames(human))]
human_DS <- human_DS[rowSums(human_DS[,2:ncol(human_DS)])>9,]
human_WT <- human[,grepl('CTRL|taxon', colnames(human))]
human_WT <- human_WT[rowSums(human_WT[,2:ncol(human_WT)])>9,]
colnames(data) <- c('taxon', 'WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')
mice_DS <- data[, grepl('DS|taxon', colnames(data))]
mice_DS <- mice_DS[rowSums(mice_DS[,2:ncol(mice_DS)])>9,]
mice_WT <- data[, grepl('WT|taxon', colnames(data))]
mice_WT <- mice_WT[rowSums(mice_WT[,2:ncol(mice_WT)])>9,]

# human
png('venn_human_DS_WT.png')
plot(euler(list('Human DS\nn = 147' = human_DS$taxon, 'Human healthy\nn = 130' = human_WT$taxon)), fills = c('#FA552F', '#599799'))
dev.off()

# humans and mice
png('venn_human_mice_DS_WT.png')
plot(euler(list('Human healthy\nn = 130' = human_WT$taxon, 'Human DS\nn = 147' = human_DS$taxon, 
                'Mice WT\nn = 354' = mice_WT$taxon, 'Mice DS\nn = 325' = mice_DS$taxon)), fills = c('#599799', '#FA552F', '#00BFC4', '#F8766D'))
dev.off()

#################### Correlate mouse with mouse abundances data and human with human, and compare them #################### 

human <- read.csv('./split_file_combined_normalized_Biagi_S.csv', stringsAsFactors = F, header=T)
human$level<-human$X <- NULL
colnames(human)<-c('taxon', paste0(rep('DOWN0', 9), 1:9), paste0(rep('DOWN', 8), 10:17), paste0(rep('CTRL0', 9), 1:9),  paste0(rep('CTRL', 7), 10:16))

data <- read.csv('split_file_combined_normalized_16S_S.csv', header=T, stringsAsFactors = F) 
data$level <- data$X <- NULL

common_species <- human$taxon[human$taxon%in%data$taxon]

library(ggplot2)
library(scales)
library(ggrepel)
library(data.table)
library(reshape2)
library(grid)
library(ggsignif)

#data <- data[data$taxon%in%common_species,]
colnames(data) <- c('taxon', 'WT1_T4', 'WT1_T0', 'WT2_T4', 'WT2_T0', 'WT3_T4', 'WT3_T0', 'WT4_T4', 'WT4_T0', 'WT5_T4', 'WT5_T0', 'WT6_T4', 'WT6_T0',
                    'DS1_T4', 'DS1_T0', 'DS2_T4', 'DS2_T0', 'DS3_T4', 'DS3_T0', 'DS4_T4', 'DS4_T0', 'DS5_T4', 'DS5_T0', 'DS6_T4', 'DS6_T0')
data_corr <- human
#data_corr <- human[human$taxon%in%common_species,]
#data_corr <- data_corr[order(match(data_corr$taxon,data$taxon)),]
x <- data.frame(matrix(data=NA, ncol = ncol(human), nrow = nrow(data)-nrow(human)))
colnames(x) <- colnames(data_corr)
data_corr <- rbind(data_corr, x)
data_corr_test <- cbind(data_corr,data[,names(data) != "taxon"])
#row.names(data_corr_test) <- data_corr_test$taxon
data_corr_test$taxon <- NULL
data_corr_test <- data_corr_test[rowMeans(data_corr_test)>19,]
# plot(1:nrow(data_corr_test), data_corr_test[,1])
# points(1:nrow(data_corr_test), data_corr_test[,2], col='red')
# axis(1, 1:nrow(data_corr_test), rownames(data_corr_test), las=2)

# Normalise:
# sums <- apply(data_corr_test, 2, sum)
# quartile <- quantile(sums)[[3]] # second quartile
# data_corr_norm <- data_corr_test
# for (i in 1:nrow(data_corr_test)){
#   data_corr_norm[i, ] <- data_corr_test[i,] * quartile / sums
# }
# data_corr_norm <- round(data_corr_norm)

WT_human <- colnames(data_corr_test)[grepl('CTRL', colnames(data_corr_test))]
mice1 <- WT_human
mice2 <- mice1 [ mice1 != mice1[1]]
WT_human_kor <- c()
WT_human_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr_test[,i], data_corr_test[,j])
      WT_human_kor <- c(WT_human_kor, k$estimate)
      WT_human_p <- c(WT_human_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

DS_human <- colnames(data_corr_test)[grepl('DOWN', colnames(data_corr_test))]
mice1 <- DS_human
mice2 <- mice1 [ mice1 != mice1[1]]
DS_human_kor <- c()
DS_human_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr_test[,i], data_corr_test[,j])
      DS_human_kor <- c(DS_human_kor, k$estimate)
      DS_human_p <- c(DS_human_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

WT_mice <- colnames(data_corr_test)[grepl('WT', colnames(data_corr_test))]
mice1<- WT_mice
mice2 <- mice1 [ mice1 != mice1[1]]
WT_mice_kor <- c()
WT_mice_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr_test[,i], data_corr_test[,j])
      WT_mice_kor <- c(WT_mice_kor, k$estimate)
      WT_mice_p <- c(WT_mice_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}

DS_mice <- colnames(data_corr_test)[grepl('DS', colnames(data_corr_test))]
mice1<- DS_mice
mice2 <- mice1 [ mice1 != mice1[1]]
DS_mice_kor <- c()
DS_mice_p <- c()
cnt<-0
for (i in mice1) {
  for (j in mice2) {
    if (i!=j){
      cat(i, ' ', j, '\n')
      k<- cor.test(data_corr_test[,i], data_corr_test[,j])
      DS_mice_kor <- c(DS_mice_kor, k$estimate)
      DS_mice_p <- c(DS_mice_p, k$p.value)
      cnt<-cnt+1 } } 
  mice2 <- mice2 [ mice2!=i ]
}


maks <- max(length(WT_human_kor), length(DS_human_kor), length(WT_mice_kor), length(DS_mice_kor))

df = data.frame('WT_human_kor' = c(WT_human_kor, rep(NA, maks-length(WT_human_kor))),
          'DS_human_kor' = c(DS_human_kor, rep(NA, maks-length(DS_human_kor))),
          'WT_mice_kor' = c(WT_mice_kor, rep(NA, maks-length(WT_mice_kor))),
          'DS_mice_kor' = c(DS_mice_kor, rep(NA, maks-length(DS_mice_kor))))


class(df)
data_plot_corr<-melt(df)

ggplot(data_plot_corr, aes(x=variable, y=value)) +
  geom_boxplot(fill = c('#7CAE00', '#C77CFF', '#00BFC4', '#F8766D'), alpha =0.8) +
  ylab('Pearson R') + xlab('Condition') +
  scale_x_discrete (labels = c('humans healthy', 'humans DS', 'mice WT', 'mice Ts65Dn')) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  geom_signif(comparisons = list(c("WT_human_kor", "DS_human_kor"),
                                 c('WT_mice_kor', 'DS_mice_kor')), 
              map_signif_level=TRUE)

ggsave('correlations_of_samples_within_genotypes_mice_and_humans.svg')

wilcox.test(c(WT_human_kor), c(DS_human_kor))
wilcox.test(c(WT_mice_kor), c(DS_mice_kor))

wilcox.test(c(WT_mice_kor, DS_mice_kor), c(WT_human_kor, DS_human_kor))

median(c(WT_human_kor, DS_human_kor))
median(c(WT_mice_kor, DS_mice_kor))
median(c(WT_human_p, DS_human_p))
median(c(WT_mice_p, DS_mice_p))
#shapiro.test(WT_human_kor)
#shapiro.test(rnorm(100, 10))

#################### Clustering of the human data alone ####################

library(vegan)
library(dendextend)

data_hcl <- human
rownames(data_hcl) <- data_hcl$taxon
data_hcl$taxon <- NULL
data_hcl <- data_hcl[rowSums(data_hcl)>19,]
data_hcl<-data.frame(t(data_hcl))
dist_matr <- vegdist(data_hcl,method="canberra") # canberra daje 3 klastry z czego 2 sa dosyc specyficzne i ward.D2 / complete poki co w miare ok
dend<-as.dendrogram(hclust(dist_matr, method='ward.D2')) ##
#labels_colors(dend)<-c(rep('#00BFC4', 8), '#F8766D', rep('#00BFC4', 3), '#F8766D', '#00BFC4', rep('#F8766D', 10))
png(file="HCL_strain_S_human.png")
#par(mar=c(5,4.1,2.1,2.1))
plot(dend)
legend("topright", c('WT', 'DS'), pch=19, col=c('#00BFC4', '#F8766D'))
dev.off()

dev.off() 

