#####################################################################################
# RNA-SEQ ANALYSIS for vascular anomalies
# Author: Beatriz Urda García, 2024
######################################################################################

# Load Libraries
library(pvclust)
library(fpc)
library(readr)
library(knitr)
library(data.table)
library(SummarizedExperiment)
library(edgeR)
#library(jaccard)
library(Rtsne)
library(tsne)
library(apcluster)
# library(Spectrum)
library(dbscan)
library(reshape2)
library(rgl)
library(ggplot2)
library(sva) # combat
library("corpcor") # svd

# Set Global Parameters
set.seed(1)
getwd()
try(setwd("COUNTS_FROM_ANE/"))
try(dev.off())

# Set Configuration
remove_controls = TRUE
new_samples = TRUE
hierarchichal_clustering = FALSE


#### Data Processing ####

### Load data
# Original counts
counts1 <- read.table('FullTable_RawCounts.txt',sep='\t'); dim(counts1)  # 58676    37
counts1[1:5, 1:5]

# New counts
counts2 <- read.table('../New_counts/FullCounts_MG-07_Illumina_TotalRNAseq.txt',sep='\t'); dim(counts2)  # 58676    18
counts2[1:5, 1:5]

# Metadata
se <- readRDS("../vascular-disease-sofia/vascular-disease-master//MG-04_Illumina_totalRNASeq/preprocess/patients_and_controls")
all_meta = as.data.frame(colData(se))
head(all_meta)

### Quality check and metadata organization
#  Transforming sample names from "VM0" followed by one or more numbers is changed to "VM" directly followed by numbers while removing initial zeros.
rownames(all_meta) = gsub("VM0*(\\d+)", "VM\\1", rownames(all_meta)) 

# Removing control samples
if(remove_controls){ # Removing 1 control in the original data set (HUVECS)
  print(dim(counts1))
  counts1 <- counts1[, !(names(counts1) %in% "HUVECS")]
  print(dim(counts1))
  print(dim(all_meta))
  meta1 = all_meta[all_meta$Case.Control == "case", ]
  print(dim(meta1))
}

# Do counts1 and meta match?
nrow(meta1) == ncol(counts1) # TRUE
N1 = nrow(meta1); N1 # 36
N1 == length(intersect(rownames(meta1),colnames(counts1))) # TRUE
N2 = ncol(counts2); N2 # 18

# Counts 1 and counts 2 have the same number of genes
nrow(counts1) == nrow(counts2) # TRUE

# Add Metadata 2
# Create an empty data frame with the specified column names
meta2 <- data.frame(matrix(ncol = ncol(meta1), nrow = N2))

# Assign the column and row names
colnames(meta2) <- colnames(meta1)
rownames(meta2) <- colnames(counts2)
head(meta2)
meta2$Case.Control = rep("case", N2)
meta2$Mutant.gene = rep("nd", N2)

# Add mutations
meta2[rownames(meta2) == "VM154", ]$Mutant.gene = "MAP3K3"
meta2[rownames(meta2) == "VM159", ]$Mutant.gene = "GNAO"
meta2[rownames(meta2) == "VM165", ]$Mutant.gene = "KRAS"
meta2[rownames(meta2) == "VM189", ]$Mutant.gene = "KRAS"
meta2[rownames(meta2) == "VM200", ]$Mutant.gene = "PIEZO1"
# meta2[rownames(meta2) == "VM204", ]$Mutant.gene = "PIK3CA;KRAS" # Finally, the variant is not present
meta1[rownames(meta1) %in% c("VM68", "VM93"), ]$Mutant.gene = "SDHD"
meta2[rownames(meta2) == "VM154", ]$Mutant.gene = "SDHD;MAP3K3"

# Create a new column 'batch' initialized with NA
meta1$batch <- NA
meta2$batch <- NA

# Assign values based on conditions in 'Process.ID' column
meta1$batch <- 1
meta2$batch <- 2

# SDHD value: TRUE, FALSE
meta1$sdhd = FALSE
meta2$sdhd = FALSE
meta1[rownames(meta1) %in% c("VM68", "VM93"), ]$sdhd = TRUE
meta2[rownames(meta2) == "VM154", ]$sdhd = TRUE

dim(counts1)
dim(meta1)
dim(counts2)
dim(meta2)

if(identical(rownames(counts1), rownames(counts2))){
  counts = cbind(counts1, counts2) 
  print(dim(counts)) # 54
}else{
  print("Genes are not sorted in the same manner.")
}

head(counts)
meta = rbind(meta1, meta2)

# Creating SummarizedExperiment object
se <- SummarizedExperiment(assays=list(counts=counts), colData=meta)
table(se$Case.Control)
table(se$batch)
table(se$batch, se$Mutant.gene)

#### RNA-seq Analysis ####

pdf(file="../Results/Vascular_plots_new_samples.pdf", width = 10, height = 10)

dge <- DGEList(counts = counts, genes = rownames(counts))
head(dge$counts)
dim(dge$counts)
write.table(dge$counts, file = '../Results/new_counts.tsv',sep= '\t', row.names = T,col.names=T, quote=F)
head(cpm(dge$counts))

# Computing counts per million
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
#head(assays(se)$logCPM)
logCPM <- cpm(dge$counts, log=TRUE, prior.count=0.5) 

## Examine sequencing depth
ord <- order(dge$samples$lib.size)
head(dge$samples$lib.size)

# Barplot patients/controls
# barplot(dge$samples$lib.size[ord]/1e+06, las = 1, main="Barplot - library size", ylab = "Million of counts", xlab = "Samples", col = c("blue","red")[(se$type[ord] == "patient") +1], border = NA)
# legend("topleft", c("patient","control"), fill = c("red","blue"), inset = 0.01)
write.table(assays(se)$logCPM, file = '../Results/new_log_CPM_counts.tsv',sep= '\t', row.names = T,col.names=T, quote=F)

# Filtering lowly expressed genes
nsamples <- length(se$Process.ID)
sample_cutoff <- 0.2
nsamples_cutoff <- sample_cutoff*nsamples

logcpm_cutoff <- 1
mask <- rowSums(logCPM <= logcpm_cutoff) >= nsamples_cutoff

dim(se)
se.filt <- se[!mask, ]; dim(se.filt)
dge.filt <- dge[!mask, ]; dim(dge.filt)
kept_genes2b <- dim(se.filt)[1]

# Between sample normalization
dge.filt.norm <- calcNormFactors(dge.filt)
assays(se.filt)$logCPM.norm <- cpm(dge.filt.norm, log = TRUE, prior.count = 3, normalized.lib.sizes = TRUE)

logCPM <- cpm(dge.filt.norm, log=TRUE, prior.count=3)   
head(logCPM)
scaled_filt_logCPM <- logCPM
write.table(file='../Results/new_scaled_expression_vascular_disease.tsv', scaled_filt_logCPM,sep='\t',row.names = T,col.names = T)
# library(WGCNA)
# dist_mat <- 1 - WGCNA::cor(scaled_filt_logCPM)
# head(dist_mat)

# Load metadata
metadata <- read.table('../COUNTS_FROM_ANE/datTraits_all.txt',sep='\t',header = T,row.names = 1)
# Filter outlier patients
#metadata <- metadata[- which(rownames(metadata) %in% c('VM024','VM040','VM099')),]
head(metadata)
rownames(metadata) <- metadata[,1]
splited_mutations <- split(rownames(metadata),metadata$Gene)
splited_mutations
splited_conditions <- split(rownames(metadata),metadata$Summary.clinic)
splited_conditions
colnames(logCPM)

# t-SNE Plotting (Colored by mutation)
library(scatterplot3d)
as.factor(metadata$Gene)
as.double(as.factor(metadata$Gene))
length(as.factor(metadata$Gene))
rownames(metadata)
length(colors)
options(repr.plot.width=20, repr.plot.height=15)
y <- Rtsne::Rtsne(t(scaled_filt_logCPM),dim=3,perplexity=10)

colores <- c('black','purple','cyan','orange','cyan')
colors <- colores[as.double(as.factor(metadata$Gene))]
# px <- scatterplot3d::scatterplot3d(y$Y[,1],y$Y[,2],y$Y[,3],angle=55,color=colors,pch=16)
# text(px$xyz.convert(y$Y[,1:3]),labels = rownames(metadata),cex= 0.7, col = "steelblue",pos = 2)
metadata
head(y$Y)


# Ward clustering
library(RColorBrewer)
options(repr.plot.width=20, repr.plot.height=10)
pal <- colorRampPalette(brewer.pal(9, "PRGn"))(100)
funky <- function(x) {as.dendrogram(hclust(x,'ward.D2'))} # Pearson Correlation
# stats::heatmap(as.matrix(dist_mat),hclustfun = funky,col=pal)
# options(device.ask.default = FALSE)

set.seed(1)

# t-SNE parameters
nsamples = nrow(meta); nsamples
max_perplex <- trunc(((nsamples - 1)/3)-1)
optimal_perplex = round(nsamples^(1/2))
#max_perplex = nsamples/3
perplex1 = round(optimal_perplex/2)
perplex3 = mean(c(optimal_perplex,max_perplex))

# t-SNEs based on logCMP normalized counts
meta$sample.name = rownames(meta)
tsne_out <- Rtsne(as.data.frame(t(assays(se.filt)$logCPM.norm)), perplexity=optimal_perplex)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$sample = meta$sample.name; 
tsne_plot = merge(tsne_plot, meta, by.x='sample', by.y='sample.name'); dim(tsne_plot)
dim(tsne_plot)
head(tsne_plot)
tsne_plot$batch = as.factor(tsne_plot$batch)

# ggplot(tsne_plot,label=batch) + geom_point(aes(x=X1,y=X2))
meta$sample.name = rownames(meta)
ggplot(tsne_plot, aes(x=X1, y=X2, color=batch, shape=Case.Control)) + geom_point() + 
  theme_classic()+
  ggtitle("t-SNE Plot based on logCMP normalized counts")
ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=batch, label=sample)) + geom_point() + 
  theme_classic()+
  ggtitle("t-SNE Plot based on logCMP normalized counts")
ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=batch, label=sample)) + geom_point() +
  geom_text(aes(label = sample), vjust = -0.5, hjust = -0.5, size = 3) +
  theme_classic() +  
  ggtitle("t-SNE Plot based on logCMP normalized counts")
# ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=Case.Control, label=sample)) + geom_point() +
#   geom_text(aes(label = sample), vjust = -0.5, hjust = -0.5, size = 3) +
#   theme_classic() +  
#   ggtitle("t-SNE Plot based on logCMP normalized counts")


# norm_clusters <- plot_hierarchichal_clustering(assays(se.filt)$logCPM.norm, 'Normalized clusters',save_clusters=TRUE, key_word="normalized", n_boot=100)

#### Differential Expression Analysis #####
mod <- model.matrix(~sdhd + batch, colData(se.filt))
mod0 <- model.matrix(~batch, colData(se.filt))
colnames(mod)
pValues <- f.pvalue(logCPM,mod,mod0)
sum(p.adjust(pValues, method="fdr") < 0.05)

# Distribution of expected p-values
par(mfrow=c(1, 1))
hist(pValues, main="Distribution of expected p-values after correcting for series", las=1)

# Mean-variance relationship
FDRcutoff <- 0.05
v <- voom(dge.filt.norm, mod, plot = TRUE) # Voom is applied to the normalized and filtered DGEList object
fit <- lmFit(v, mod)
fit<- eBayes(fit)
res <- decideTests(fit, p.value = FDRcutoff)
summary(res)

# Add gene metadata (gene symbol) and fetch the table of results. Print first 15 sDEgenes
rowRanges(se.filt)
#genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se.filt))), symbol = rowData(se.filt)[, 1], stringsAsFactors = FALSE)
genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt, 15)
write.table(tt,file="../Results/DEG_sdhd_vs_other.txt",sep="\t",col.names=T,row.names=F, quote=FALSE)

sDEGs = tt[tt$adj.P.Val < 0.05, ]
write.table(sDEGs,file="../Results/sDEG_sdhd_vs_other.txt",sep="\t",col.names=T,row.names=F, quote=FALSE)

# Total number of DE genes (considering p-value)
sum(tt$P.Value < 0.05)

# When considering adjusted p-value, the number os significantly DEgenes is:
sum(tt$adj.P.Val < 0.05)

# Volcano plot
par(xpd=T, mar=par()$mar+c(0,0,0,5))
with(tt, plot(logFC, -log10(adj.P.Val), pch = 20, main = "", xlim=c(-6,7), ylim=c(0,15))) 


with(subset(tt, logFC > 1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="green"))
with(subset(tt, logFC < -1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="orange"))
with(subset(tt, -log10(adj.P.Val) < 1.3), points(logFC, -log10(adj.P.Val), pch = 20, col="grey"))

legend("bottomright", cex = .75, inset = c(-0.32,0.75), xjust = 2, yjust = 1,pch = 20,c("sDE & overexpressed", "sDE & underexpressed", "non-sDE", "sDE"), col = c("green", "orange", "grey", "black"),bg="white", bty= "n")
with(subset(tt, -log10(adj.P.Val)>1.3 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs = symbol, cex=.7, col= "#35445b"))
par(xpd=F)
abline(h= 1.3, col = "blue", lty= 2, lwd = 1)
abline(v= c(-1,1), col = "blue", lty= 2, lwd = 1)

# Select genes with a |logFC| higher than 1
sum((tt$adj.P.Val < 0.05) & (abs(tt$logFC)) >= 1)

#### Batch Effect Correction ####

### Combat 
batch <- se.filt$batch
# without protecting biological variables of interest
combatexp <- ComBat(assays(se.filt)$logCPM.norm, batch, mod=NULL) # matrix with the batch effect corrected
class(combatexp)
dim(combatexp)
# protecting the effect of mutations
mod <- model.matrix(~as.factor(Mutant.gene), colData(se.filt))
combatexp <- ComBat(assays(se.filt)$logCPM.norm, batch, mod=mod, par.prior = TRUE) # matrix with the batch effect corrected
write.table(file='../Results/combat_expression_vascular_disease.csv', as.data.frame(combatexp),sep='\t',row.names = T, col.names = T)

# t-SNEs based on Combat corrected logCMP normalized counts
tsne_out <- Rtsne(as.data.frame(t(combatexp)), perplexity=optimal_perplex)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$sample = meta$sample.name; 
tsne_plot = merge(tsne_plot, meta, by.x='sample', by.y='sample.name'); dim(tsne_plot)
tsne_plot$batch = as.factor(tsne_plot$batch)
dim(tsne_plot)

ggplot(tsne_plot, aes(x=X1, y=X2, color=batch, shape=Case.Control)) + geom_point() + 
  theme_classic()+
  ggtitle("t-SNE Plot based on Combat corrected counts")
ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=as.character(batch), label=sample)) + geom_point()+
  theme_classic()+
  ggtitle("t-SNE Plot based on Combat corrected counts")

ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=as.character(batch), label=sample)) + geom_point() +
  geom_text(aes(label = sample), vjust = -0.5, hjust = -0.5, size = 3) +
  theme_classic()+
  ggtitle("t-SNE Plot based on Combat corrected counts")


### SVD - it removes unknown batch effect, it behaves worse
s <- fast.svd(t(scale(t(assays(se.filt)$logCPM.norm), center = TRUE, scale = TRUE)))
pcSds <- s$d
pcSds[1] <- 0
svdexp <- s$u %*% diag(pcSds) %*% t(s$v)
colnames(svdexp) <- colnames(se.filt)
class(svdexp)
dim(svdexp)

d <- as.dist(1 - cor(svdexp, method = "spearman"))
sampleClustering <- hclust(d)

# t-SNEs based on SVD corrected logCMP normalized counts
tsne_out <- Rtsne(as.data.frame(t(svdexp)), perplexity=optimal_perplex)
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$sample = meta$sample.name; 
tsne_plot = merge(tsne_plot, meta, by.x='sample', by.y='sample.name'); dim(tsne_plot)
tsne_plot$batch = as.factor(tsne_plot$batch)
dim(tsne_plot)

ggplot(tsne_plot, aes(x=X1, y=X2, color=batch, shape=Case.Control)) + geom_point() + 
  theme_classic()+
  ggtitle("t-SNE Plot based on SVD corrected counts")
ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=as.character(batch), label=sample)) + geom_point()+
  theme_classic()+
  ggtitle("t-SNE Plot based on SVD corrected counts")
ggplot(tsne_plot, aes(x=X1, y=X2, color=Mutant.gene, shape=as.character(batch), label=sample)) + geom_point() +
  geom_text(aes(label = sample), vjust = -0.5, hjust = -0.5, size = 3) +
  theme_classic()+
  ggtitle("t-SNE Plot based on SVD corrected counts")

# svd_clusters <- plot_hierarchichal_clustering(svdexp, 'after SVD batch effect removal',save_clusters=TRUE, key_word="svd", n_boot=100)
# combat_clusters <- plot_hierarchichal_clustering(combatexp, 'after Combat batch effect removal',save_clusters=TRUE, key_word="combat", n_boot=100)

if(hierarchichal_clustering){
  library(pvclust)
  library(dendextend)
  library(ggdendro)
  result <- pvclust(combatexp, method.hclust="ward.D2", method.dist="euclidean", nboot=100)
  
  # Print the result
  print(result)
  
  # Plot the dendrogram with p-values
  plot(result)
  pvrect(result, alpha=0.95) # Add rectangles around clusters with approximately unbiased (AU) p-values greater than 0.95
  
  # Obtener el dendrograma
  dend <- as.dendrogram(result$hclust)
  
  # Añadir los metadatos al dendrograma
  dend <- dend %>% 
    set("labels", meta$Mutant.gene[match(labels(dend), meta$sample.name)])
  
  # Convertir el dendrograma a un data frame para ggplot2
  dend_data <- dendro_data(dend, type = "rectangle")
  
  # Añadir la información de los metadatos al data frame
  dend_data$labels <- merge(dend_data$labels, meta, by.x = "label", by.y = "sample.name")
  
  # Crear el gráfico con ggplot2
  p <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = dend_data$labels, aes(x = x, y = y, label = label, color = Mutant.gene), size = 3, hjust = 1) +
    scale_color_manual(values = c("SDHD" = "red", "PIK3CA" = "blue", "nd" = "green")) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # Mostrar el gráfico
  print(p)
  
}


dev.off()



