library(ggplot2)
library(gridExtra)
library(DESeq2)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(RColorBrewer)
library(factoextra)
library(Seurat)
library(readxl)
library(ggpubr)


### All over the script: 
### spi >>> Late stationary phase
### ana >>> Anaerobic shock
### salt >>> Salt shock


      
setwd("/address/to/working/directory/")
dir.create("output")
            
## importing col Data
colData <- read.csv("./lib_changed_corr_fa.csv", row.names = 1)
rownames(colData) <- as.character(colData$library)
## importing row data
anno <- read.table("./index/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.ASM21085v2.45.gtf", sep = "\t")
anno <- anno[anno$V3 == "gene", ]
            
biotype <- gsub(";", "", gsub(".*; gene_biotype ", "", anno$V9))
gene_name <- gsub("; .*", "", gsub("gene_id.*", "",gsub(".*; gene_name ", "", anno$V9)))
gene_id <- gsub("; .*", "",gsub("gene_id ", "", anno$V9))
            
rowData <- data.frame(gene_name, biotype)
rownames(rowData) <- gene_id
            
## Importing the count table
temp <- list.files("./unique_count", pattern = "unique.txt")
col <- gsub("_S.*", "", gsub("Aligned.out.unique.txt", "", temp))
        
filenames <- list.files(path = "./unique_count", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,6]})
unique_table <- sapply(datalist, cbind)
colnames(unique_table) <- col
rownames(unique_table) <- rownames(rowData)
          
## limiting data to single and ten cells library
colData <- colData[colData$n.cells == "single" | colData$n.cells == "ten",]
colData <- colData[colnames(unique_table),][-c(90,113,51,54,102), ]
unique_table <- unique_table[,rownames(colData)]
          
## Adding the number of detected genes and library size
colData$n.gene <- colSums(unique_table > 5)
colData$lib.size <- colSums(unique_table)
          

            
## Fig2b
p1 <- ggplot(data = colData[colData$n.cells == "single", ], aes(x = condition, y=n.gene, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  ggtitle("Single Cell")
            
            
p2 <- ggplot(data = colData[colData$n.cells == "ten", ], aes(x = condition, y=n.gene, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") 
Fig2b <- grid.arrange(p1, p2, ncol=1)
ggsave("./output/Fig2b.png", Fig2b, width = 8, height = 6)
            
            
## sFig2b
p1 <- ggplot(data = colData[colData$n.cells == "single", ], aes(x = lib.size, y=n.gene, color = condition)) +
  geom_jitter() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Number of detected genes") +
  xlab("Library Size") +
  ggtitle("single cell")
              
p2 <- ggplot(data = colData[colData$n.cells == "ten", ], aes(x = lib.size, y=n.gene, color = condition)) +
  geom_jitter() +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Number of detected genes") +
  xlab("Library Size") +
  ggtitle("10 pooled cells")
        
sFig2b <- grid.arrange(p1, p2, ncol=2)
ggsave("./output/sFig2b.png", sFig2b, width = 10, height = 5)
          
          
## sFig2a
p1 <- ggplot(data = colData[colData$n.cells == "single", ], aes(x = condition, y=lib.size, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  ggtitle("Single Cell")
            
            
p2 <- ggplot(data = colData[colData$n.cells == "ten", ], aes(x = condition, y=lib.size, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  ggtitle("10 pooled cells")
            
sFig2a <- grid.arrange(p1, p2, ncol=2)
ggsave("./output/sFig2a.png", sFig2a, width = 12, height = 4)
          
          
          
single <- colData[colData$n.cells == "single",]
ten <- colData[colData$n.cells == "ten",]
          
          
single_table <- unique_table[,rownames(single)]
ten_table <- unique_table[,rownames(ten)]
          
          
### Normalizing data
sfSingle <- estimateSizeFactorsForMatrix(single_table)
sfTen <- estimateSizeFactorsForMatrix(ten_table)
          
nSingle <- t(t(single_table) / sfSingle)
nTen <- t(t(ten_table) / sfTen)
          
##########################################################################
## Selecting the top variable genes and log transformation
topVarGenesSingle <- head(order(rowVars(nSingle), decreasing = TRUE),300)
nSingle_fin <- log(nSingle[topVarGenesSingle,] + 1)
          
topVarGenesTen <- head(order(rowVars(nTen), decreasing = TRUE),300)
nTen_fin <- log(nTen[topVarGenesTen,] + 1)
          
          
### Performing PCA (ten cells)
res.pca <- prcomp(t(nTen_fin), scale = FALSE)
pca_ten <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = ten$condition), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  ggtitle("Ten cells")

## sFig6b            
lib_pca_ten <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = ten$lib.size), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  scale_color_gradientn(colours = brewer.pal(9, "YlOrRd")) +
  theme_classic() +
  theme() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  ggtitle("lib_size (Ten cells) ")
          
          
n.gene_pca_ten <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = ten$n.gene), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  scale_color_gradientn(colours = brewer.pal(9, "YlOrRd")) +
  theme_classic() +
  theme() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  ggtitle("n.gene (Ten cells)")
        
sFig6b <- grid.arrange(lib_pca_ten, n.gene_pca_ten, ncol=2)
ggsave("./output/sFig6b.png", sFig6b, width = 10)
          
## sFig6c_left
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 20))
ggsave("./output/sFig6c_left.png")
 
## sFig6b_biplot
var <- get_pca_var(res.pca)
contrib <- fviz_contrib(res.pca, choice = "var", axes = 1:2)$data
contrib <- rownames(contrib[order((contrib$contrib), decreasing = TRUE), ])[1:15]


biplot_ten <- fviz_pca_biplot(res.pca, repel = TRUE, geom = c("point"), geom.var = c("arrow"), addEllipses = TRUE, 
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969",  # Individuals color
                select.var = list(name = contrib)
)
biplot_ten <- biplot_ten +  
  geom_jitter(data  = data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = ten$condition), aes(x = PC1, y = PC2, color = condition)) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("PC1") +
  ylab("PC2")
  
ggsave("./output/sFig6b_biplot.png", height = 4, width = 3.5)


### Performing PCA (single cell)
res.pca <- prcomp(t(nSingle_fin), scale = FALSE)
pca_single <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = single$condition), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  ggtitle("Single cell")

##sFig6a
lib_pca_single <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = single$lib.size), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  scale_color_gradientn(colours = brewer.pal(9, "YlOrRd")) +
  theme_classic() +
  theme() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  ggtitle("lib_size (single cells) ")
          
          
n.gene_pca_single <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = single$n.gene), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  scale_color_gradientn(colours = brewer.pal(9, "YlOrRd")) +
  theme_classic() +
  theme() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  ggtitle("n.gene (single cells)")
          
          
sFig6a <- grid.arrange(lib_pca_single, n.gene_pca_single, ncol=2)
ggsave("./output/sFig6a.png", sFig5a, width = 10)
 

## Fig2a         
Fig3a <- grid.arrange(pca_single, pca_ten, ncol=1)
ggsave("./output/Fig3a.png", Fig3a, height = 12)
          
## sFig6c_right
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 20))
ggsave("./output/sFig6c_right.png")

## sFig6a_biplot
var <- get_pca_var(res.pca)
contrib <- fviz_contrib(res.pca, choice = "var", axes = 1:2)$data
contrib <- rownames(contrib[order((contrib$contrib), decreasing = TRUE), ])[1:15]


biplot_single <- fviz_pca_biplot(res.pca, repel = TRUE, geom = c("point"), geom.var = c("arrow"), addEllipses = TRUE, 
                              col.var = "#2E9FDF", # Variables color
                              col.ind = "#696969",  # Individuals color
                              select.var = list(name = contrib)
)
biplot_single <- biplot_single +  
  geom_jitter(data  = data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = single$condition), aes(x = PC1, y = PC2, color = condition)) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("PC1") +
  ylab("PC2")

ggsave("./output/sFig6a_biplot.png", height = 4, width = 3)
#############################################################################
################ Correlation between our data and Kroger data ###############
#############################################################################
### reading the table from Kroger data
tpm <- read_excel("./Kroger/mmc3.xlsx")
colnames(tpm) <- tpm[2,]
tpm <- tpm[-c(1,2), ]
tpm <- as.data.frame(tpm)
rownames(tpm) <- tpm[,1]

## sFig5_salt_ten
salt_ten_pooled <- data.frame(mean = rowMeans(nTen[ ,rownames(ten[ten$condition == "salt",])]), gene = rownames(nTen))
rownames(salt_ten_pooled) <- gsub("1344_", "", rownames(salt_ten_pooled))
salt_ten_pooled <- salt_ten_pooled[rownames(tpm),]
salt_ten_pooled <- data.frame(salt_ten_pooled, bulk = as.numeric(tpm[,19]))

ggplot(salt_ten_pooled, aes(x = bulk + 1, y = mean + 1)) +
  geom_point(color = "#00BA38") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean expression in our data (log)") +
  xlab("tpm in Hinton data (log)") +
  stat_cor(method = "spearman") +
  theme_classic()

ggsave("./output/sFig5_salt_ten.png")

## sFig5_ana_ten 
ana_ten_pooled <- data.frame(mean = rowMeans(nTen[ ,rownames(ten[ten$condition == "ana",])]), gene = rownames(nTen))
rownames(ana_ten_pooled) <- gsub("1344_", "", rownames(ana_ten_pooled))
ana_ten_pooled <- ana_ten_pooled[rownames(tpm),]
ana_ten_pooled <- data.frame(ana_ten_pooled, bulk = as.numeric(tpm[,22]))


ggplot(ana_ten_pooled, aes(x = bulk + 1, y = mean + 1)) +
  geom_point(color = "#F8766D") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean expression in our data (log)") +
  xlab("tpm in Hinton data (log)") +
  stat_cor(method = "spearman") +
  theme_classic()

ggsave("./output/sFig5_ana_ten.png")

## sFig5_spi_ten
spi_ten_pooled <- data.frame(mean = rowMeans(nTen[ ,rownames(ten[ten$condition == "spi",])]), gene = rownames(nTen))
rownames(spi_ten_pooled) <- gsub("1344_", "", rownames(spi_ten_pooled))
spi_ten_pooled <- spi_ten_pooled[rownames(tpm),]
spi_ten_pooled <- data.frame(spi_ten_pooled, bulk = as.numeric(tpm[,14]))


ggplot(spi_ten_pooled, aes(x = bulk + 1, y = mean + 1)) +
  geom_point(color = "#619cff") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean expression in our data (log)") +
  xlab("tpm in Hinton data (log)") +
  stat_cor(method = "spearman") +
  theme_classic()

ggsave("./output/sFig5_spi_ten.png")


## sFig5_salt_single
salt_single_pooled <- data.frame(mean = rowMeans(nSingle[ ,rownames(single[single$condition == "salt",])]), gene = rownames(nSingle))
rownames(salt_single_pooled) <- gsub("1344_", "", rownames(salt_single_pooled))
salt_single_pooled <- salt_single_pooled[rownames(tpm),]
salt_single_pooled <- data.frame(salt_single_pooled, bulk = as.numeric(tpm[,19]))


ggplot(salt_single_pooled, aes(x = bulk + 1, y = mean + 1)) +
  geom_point(color = "#00BA38") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean expression in our data (log)") +
  xlab("tpm in Hinton data (log)") +
  stat_cor(method = "spearman") +
  theme_classic()

ggsave("./output/sFig5_salt_single.png")

## sFig5_ana_single
ana_single_pooled <- data.frame(mean = rowMeans(nSingle[ ,rownames(single[single$condition == "ana",])]), gene = rownames(nSingle))
rownames(ana_single_pooled) <- gsub("1344_", "", rownames(ana_single_pooled))
ana_single_pooled <- ana_single_pooled[rownames(tpm),]
ana_single_pooled <- data.frame(ana_single_pooled, bulk = as.numeric(tpm[,22]))


ggplot(ana_single_pooled, aes(x = bulk + 1, y = mean + 1)) +
  geom_point(color = "#F8766D") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean expression in our data (log)") +
  xlab("tpm in Hinton data (log)") +
  stat_cor(method = "spearman") +
  theme_classic()

ggsave("./output/sFig5_ana_single.png")

## sFig5_spi_single
spi_single_pooled <- data.frame(mean = rowMeans(nSingle[ ,rownames(single[single$condition == "spi",])]), gene = rownames(nSingle))
rownames(spi_single_pooled) <- gsub("1344_", "", rownames(spi_single_pooled))
spi_single_pooled <- spi_single_pooled[rownames(tpm),]
spi_single_pooled <- data.frame(spi_single_pooled, bulk = as.numeric(tpm[,14]))


ggplot(spi_single_pooled, aes(x = bulk + 1, y = mean + 1)) +
  geom_point(color = "#619cff") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean expression in our data (log)") +
  xlab("tpm in Hinton data (log)") +
  stat_cor(method = "spearman") +
  theme_classic()

ggsave("./output/sFig5_spi_single.png")
#########################################################################
############## Correlation between single and ten data ##################
#########################################################################
salt_pooled <- data.frame(mean_ten = rowMeans(nTen[ ,rownames(ten[ten$condition == "salt",])]), mean_single = rowMeans(nSingle[ ,rownames(single[single$condition == "salt",])]))
ana_pooled <- data.frame(mean_ten = rowMeans(nTen[ ,rownames(ten[ten$condition == "ana",])]), mean_single = rowMeans(nSingle[ ,rownames(single[single$condition == "ana",])]))
spi_pooled <- data.frame(mean_ten = rowMeans(nTen[ ,rownames(ten[ten$condition == "spi",])]), mean_single = rowMeans(nSingle[ ,rownames(single[single$condition == "spi",])]))

## sFig5a_salt
ggplot(salt_pooled, aes(x = mean_single + 1, y = mean_ten + 1)) +
  geom_point(color = "#00BA38") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean ten (log)") +
  xlab("mean single (log)") +
  stat_cor(method = "pearson") +
  theme_classic()

ggsave("./output/sFig5a_salt.png")

## sFig5a_ana
ggplot(ana_pooled, aes(x = mean_single + 1, y = mean_ten + 1)) +
  geom_point(color = "#F8766D") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean ten (log)") +
  xlab("mean single (log)") +
  stat_cor(method = "pearson") +
  theme_classic()

ggsave("./output/sFig5a_ana.png")

## sFig5a_spi
ggplot(spi_pooled, aes(x = mean_single + 1, y = mean_ten + 1)) +
  geom_point(color = "#619cff") +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  ylab("mean ten (log)") +
  xlab("mean single (log)") +
  stat_cor(method = "pearson") +
  theme_classic()

ggsave("./output/sFig5a_spi.png")

########################################################
########################################################
########################################################
## sFig4a
spi_ten <- nTen[,rownames(ten[ten$condition == "spi",])]
cv_spi <- data.frame(cv = rowSds(spi_ten)/rowMeans(spi_ten) , mean = rowMeans(spi_ten))
plot_spi <- ggplot(cv_spi, aes(x = mean , y = cv)) +
  geom_point(color = "#619cff") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous() +
  theme_classic()
  
salt_ten <- nTen[,rownames(ten[ten$condition == "salt",])]
cv_salt <- data.frame(cv = rowSds(salt_ten)/rowMeans(salt_ten) , mean = rowMeans(salt_ten))
plot_salt <- ggplot(cv_salt, aes(x = mean , y = cv)) +
  geom_point(color = "#00BA38") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous() +
  theme_classic()
  
ana_ten <- nTen[,rownames(ten[ten$condition == "ana",])]
cv_ana <- data.frame(cv = rowSds(ana_ten)/rowMeans(ana_ten) , mean = rowMeans(ana_ten))
plot_ana <- ggplot(cv_ana, aes(x = mean , y = cv)) +
  geom_point(color = "#F8766D") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous() +
  theme_classic()
    
sFig4a <- grid.arrange(plot_salt, plot_ana, plot_spi, ncol = 3)
ggsave("./output/sFig4a.png", sFig4a, width = 20, height = 6)
    
    
## sFig4b
salt_single <- nSingle[,rownames(single[single$condition == "salt",])]
cv_salt <- data.frame(cv = rowSds(salt_single)/rowMeans(salt_single) , mean = rowMeans(salt_single))
plot_salt <- ggplot(cv_salt, aes(x = mean , y = cv)) +
  geom_point(color = "#00BA38") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous() +
  theme_classic()
    
ana_single <- nSingle[,rownames(single[single$condition == "ana",])]
cv_ana <- data.frame(cv = rowSds(ana_single)/rowMeans(ana_single) , mean = rowMeans(ana_single))
plot_ana <- ggplot(cv_ana, aes(x = mean , y = cv)) +
  geom_point(color = "#F8766D") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous() +
  theme_classic()
  
spi_single <- nSingle[,rownames(single[single$condition == "spi",])]
cv_spi <- data.frame(cv = rowSds(spi_single)/rowMeans(spi_single) , mean = rowMeans(spi_single))
plot_spi <- ggplot(cv_spi, aes(x = mean , y = cv)) +
  geom_point(color = "#619cff") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous() +
  theme_classic()

sFig4b <- grid.arrange(plot_salt, plot_ana, plot_spi, ncol = 3)
ggsave("./output/sFig4b.png", sFig4b, width = 20, height = 6)
    
##############################################
##############################################
##############################################
## Fig1a
### Here the non_unique counts are used to calculate the percentage of rRNA, tRNA, etc.
filenames <- list.files(path = "./non_unique_count", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,6]})
non_unique_table <- sapply(datalist, cbind)
colnames(non_unique_table) <- col
rownames(non_unique_table) <- rownames(rowData)
    
    
non_unique_table_single <- non_unique_table[,rownames(single)]
non_unique_table_ten <- non_unique_table[,rownames(ten)]
  
    
## single cells
ratio_single <- data.frame(
  rRNA = colSums(non_unique_table_single[grep("rRNA", biotype), ]),
  tRNA = colSums(non_unique_table_single[grep("tRNA", biotype), ]),
  ncRNA = colSums(non_unique_table_single[grep("ncRNA", biotype), ]),
  sRNA = colSums(non_unique_table_single[grep("sRNA", biotype), ]),
  protein_coding = colSums(non_unique_table_single[grep("protein_coding", biotype), ]),
  pseudogene = colSums(non_unique_table_single[grep("pseudogene", biotype), ]),
  antisense = colSums(non_unique_table_single[grep("antisense", biotype), ]),
  ribozyme = non_unique_table_single[grep("ribozyme", biotype), ],
  antitoxin = colSums(non_unique_table_single[grep("antitoxin", biotype), ])
  )
    
write.csv(ratio_single, file = "./output/ratio_single.csv")
    
    
percentage_single <- data.frame(
  rRNA = (colSums(non_unique_table_single[grep("rRNA", biotype), ]) / rowSums(ratio_single)) * 100,
  tRNA = (colSums(non_unique_table_single[grep("tRNA", biotype), ]) / rowSums(ratio_single)) * 100,
  ncRNA = (colSums(non_unique_table_single[grep("ncRNA", biotype), ]) / rowSums(ratio_single)) * 100,
  sRNA = (colSums(non_unique_table_single[grep("sRNA", biotype), ])  / rowSums(ratio_single)) * 100,
  protein_coding = (colSums(non_unique_table_single[grep("protein_coding", biotype), ])  / rowSums(ratio_single)) * 100,
  pseudogene = (colSums(non_unique_table_single[grep("pseudogene", biotype), ])  / rowSums(ratio_single)) * 100,
  antisense = (colSums(non_unique_table_single[grep("antisense", biotype), ])  / rowSums(ratio_single)) * 100,
  ribozyme = (non_unique_table_single[grep("ribozyme", biotype), ]  / rowSums(ratio_single)) * 100,
  antitoxin = (colSums(non_unique_table_single[grep("antitoxin", biotype), ]) / rowSums(ratio_single)) * 100
  )
    
write.csv(percentage_single, file = "./output/percentage_single.csv")

## ten cells
ratio_ten <- data.frame(
  rRNA = colSums(non_unique_table_ten[grep("rRNA", biotype), ]),
  tRNA = colSums(non_unique_table_ten[grep("tRNA", biotype), ]),
  ncRNA = colSums(non_unique_table_ten[grep("ncRNA", biotype), ]),
  sRNA = colSums(non_unique_table_ten[grep("sRNA", biotype), ]),
  protein_coding = colSums(non_unique_table_ten[grep("protein_coding", biotype), ]),
  pseudogene = colSums(non_unique_table_ten[grep("pseudogene", biotype), ]),
  antisense = colSums(non_unique_table_ten[grep("antisense", biotype), ]),
  ribozyme = non_unique_table_ten[grep("ribozyme", biotype), ],
  antitoxin = colSums(non_unique_table_ten[grep("antitoxin", biotype), ])
  )
    
write.csv(ratio_ten, file = "./output/ratio_ten.csv")

    
    
percentage_ten <- data.frame(
  rRNA = (colSums(non_unique_table_ten[grep("rRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
  tRNA = (colSums(non_unique_table_ten[grep("tRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
  ncRNA = (colSums(non_unique_table_ten[grep("ncRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
  sRNA = (colSums(non_unique_table_ten[grep("sRNA", biotype), ])  / rowSums(ratio_ten)) * 100,
  protein_coding = (colSums(non_unique_table_ten[grep("protein_coding", biotype), ])  / rowSums(ratio_ten)) * 100,
  pseudogene = (colSums(non_unique_table_ten[grep("pseudogene", biotype), ])  / rowSums(ratio_ten)) * 100,
  antisense = (colSums(non_unique_table_ten[grep("antisense", biotype), ])  / rowSums(ratio_ten)) * 100,
  ribozyme = (non_unique_table_ten[grep("ribozyme", biotype), ]  / rowSums(ratio_ten)) * 100,
  antitoxin = (colSums(non_unique_table_ten[grep("antitoxin", biotype), ]) / rowSums(ratio_ten)) * 100
  )
    
write.csv(percentage_ten, file = "./output/percentage_ten.csv")
rm(list = setdiff(ls(), c("unique_table","ten_table", "ten", "nTen", "single_table", "single", "nSingle", "colData", "rowData", "norm")))
##################################################
############# Saturation Analysis ################
##################################################
## Ten cells
lib <- seq(0, max(colSums(ten_table)), by = 1000)
gene <- data.frame(matrix(0, nrow = 57, ncol = length(lib)))
rownames(gene) <- colnames(ten_table)
        
for (i in 1:length(lib)) {
  matrix <- SampleUMI(as.matrix(ten_table), max.umi = lib[i], upsample = FALSE)
  gene[,i] <- colSums(matrix > 5)
  }
        
gene <- t(rbind(lib, gene))
colnames(gene)[1] <- "lib"
    
for (i in 1:57) {
  size <- ten$lib.size[i]
  size <- gene[,1] > size
  size <- which(size)
  gene[size, i+1] <- NA
  }
        
        
salt_gene <- gene[,rownames(ten[ten$condition == "salt", ])]
pos <- which(rowSums(is.na(salt_gene)) < (length(colnames(salt_gene)) - 5))
salt_gene <- salt_gene[pos,]
df_salt <- data.frame(salt_mean = rowMeans(salt_gene, na.rm = TRUE), salt_std = rowSds(salt_gene, na.rm = TRUE), lib = lib[pos])
      
ana_gene <- gene[,rownames(ten[ten$condition == "ana", ])]
pos <- which(rowSums(is.na(ana_gene)) < (length(colnames(ana_gene)) - 5))
ana_gene <- ana_gene[pos,]
df_ana <- data.frame(ana_mean = rowMeans(ana_gene, na.rm = TRUE), ana_std = rowSds(ana_gene, na.rm = TRUE), lib = lib[pos])
        
spi_gene <- gene[,rownames(ten[ten$condition == "spi", ])]
pos <- which(rowSums(is.na(spi_gene)) < (length(colnames(spi_gene)) - 5))
spi_gene <- spi_gene[pos,]
df_spi <- data.frame(spi_mean = rowMeans(spi_gene, na.rm = TRUE), spi_std = rowSds(spi_gene, na.rm = TRUE), lib = lib[pos])
        
## Fig2c_ten        
ggplot(df_spi) + 
  geom_smooth(data = df_spi + 1, aes(lib, spi_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#619cff", n=1000) +
  geom_smooth(data = df_ana + 1, aes(lib, ana_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#F8766D", n=1000) +
  geom_smooth(data = df_salt + 1, aes(lib, salt_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#00BA38", n=1000) +
  xlab("mean library size") +
  ylab("mean number of detected genes") +
  xlim(0,1000000) + 
  ylim(0, NA) +
  theme_classic()
    
ggsave("./output/Fig2c_ten.png", width = 6, height = 3)
   
        
# Single cell
lib <- seq(0, max(colSums(single_table)), by = 1000)
gene <- data.frame(matrix(0, nrow = 69, ncol = length(lib)))
rownames(gene) <- colnames(single_table)
        
for (i in 1:length(lib)) {
  matrix <- SampleUMI(as.matrix(single_table), max.umi = lib[i], upsample = FALSE)
  gene[,i] <- colSums(matrix > 5)
  }
        
gene <- t(rbind(lib, gene))
colnames(gene)[1] <- "lib"
        
for (i in 1:69) {
  size <- single$lib.size[i]
  size <- gene[,1] > size
  size <- which(size)
  gene[size, i+1] <- NA
  }
        
        
salt_gene <- gene[,rownames(single[single$condition == "salt", ])]
pos <- which(rowSums(is.na(salt_gene)) < (length(colnames(salt_gene)) - 5))
salt_gene <- salt_gene[pos,]
df_salt <- data.frame(salt_mean = rowMeans(salt_gene, na.rm = TRUE), salt_std = rowSds(salt_gene, na.rm = TRUE), lib = lib[pos])
        
ana_gene <- gene[,rownames(single[single$condition == "ana", ])]
pos <- which(rowSums(is.na(ana_gene)) < (length(colnames(ana_gene)) - 5))
ana_gene <- ana_gene[pos,]
df_ana <- data.frame(ana_mean = rowMeans(ana_gene, na.rm = TRUE), ana_std = rowSds(ana_gene, na.rm = TRUE), lib = lib[pos])
        
spi_gene <- gene[,rownames(single[single$condition == "spi", ])]
pos <- which(rowSums(is.na(spi_gene)) < (length(colnames(spi_gene)) - 5))
spi_gene <- spi_gene[pos,]
df_spi <- data.frame(spi_mean = rowMeans(spi_gene, na.rm = TRUE), spi_std = rowSds(spi_gene, na.rm = TRUE), lib = lib[pos])
    
ggplot(df_spi) + 
  geom_smooth(data = df_spi + 1, aes(lib, spi_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#619cff", n=1000, fullrange=TRUE, linetype = "dashed") +
  geom_smooth(data = df_spi + 1, aes(lib, spi_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#619cff", n=1000) +
  geom_smooth(data = df_ana + 1, aes(lib, ana_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#F8766D", n=1000) +
  geom_smooth(data = df_salt + 1, aes(lib, salt_mean),se = FALSE, method = "gam", formula = y ~ log(x), size = 1.5, colour = "#00BA38", n=1000) +
  xlab("mean library size") +
  ylab("mean number of detected genes") +
  xlim(0,1000000) + 
  ylim(0, NA) +
  theme_classic()

## Fig2c_single
ggsave("./output/Fig2_single.png", width = 6, height = 3)   

#################################################################################
################## using DEseq2 - one vs other - single cell ####################
#################################################################################
single$SALT <- factor(as.character(gsub("spi", "other",gsub("ana", "other", single$condition))))
single$ANA <- factor(as.character(gsub("spi", "other",gsub("salt", "other", single$condition))))
single$SPI <- factor(as.character(gsub("salt", "other",gsub("ana", "other", single$condition))))
    
#SALT
cts <- cbind(single_table[,rownames(single[single$condition == "salt",])], single_table[,rownames(single[single$condition == "ana",])], single_table[,rownames(single[single$condition == "spi",])])
cts <- cts[rowSums(cts) > 0, ]
cts <- cts[rowSums(cts > 0) > 5, ]
  
col_data <- single[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design= ~SALT)
dds <- DESeq(dds)
resultsNames(dds)
res_salt_vs_other <- results(dds, name = "SALT_salt_vs_other")
res_salt_vs_other <- data.frame(res_salt_vs_other)
res_salt_vs_other$gene_name <- rowData[rownames(res_salt_vs_other), ][,1] 
write.csv(res_salt_vs_other, file = "./output/res_salt_vs_other.csv")  
    
#ANA
cts <- cbind(single_table[,rownames(single[single$condition == "ana",])], single_table[,rownames(single[single$condition == "salt",])], single_table[,rownames(single[single$condition == "spi",])])
cts <- cts[rowSums(cts) > 0, ]
cts <- cts[rowSums(cts > 0) > 5, ]
    
col_data <- single[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design= ~ANA)
dds <- DESeq(dds)
resultsNames(dds)
res_ana_vs_other <- results(dds, name = "ANA_other_vs_ana")
res_ana_vs_other <- data.frame(res_ana_vs_other)
res_ana_vs_other$gene_name <- rowData[rownames(res_ana_vs_other), ][,1] 
write.csv(res_ana_vs_other, file = "./output/res_ana_vs_other.csv")  
    

#SPI
cts <- cbind(single_table[,rownames(single[single$condition == "spi",])], single_table[,rownames(single[single$condition == "salt",])], single_table[,rownames(single[single$condition == "ana",])])
cts <- cts[rowSums(cts) > 0, ]
cts <- cts[rowSums(cts > 0) > 5, ]

col_data <- single[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design= ~SPI)
dds <- DESeq(dds)
resultsNames(dds)
res_spi_vs_other <- results(dds, name = "SPI_spi_vs_other")
res_spi_vs_other <- data.frame(res_spi_vs_other)
res_spi_vs_other$gene_name <- rowData[rownames(res_spi_vs_other), ][,1] 
write.csv(res_spi_vs_other, file = "./output/res_spi_vs_other.csv")
    
DE_salt <- res_salt_vs_other[complete.cases(res_salt_vs_other), ]
DE_salt <- DE_salt[DE_salt$padj < 0.01,]
DE_salt <- DE_salt[order((DE_salt$log2FoldChange), decreasing = TRUE), ]
DE_salt <- DE_salt[DE_salt$log2FoldChange > 0, ]
  
DE_ana <- res_ana_vs_other[complete.cases(res_ana_vs_other), ]
DE_ana <- DE_ana[DE_ana$padj < 0.01,]
DE_ana <- DE_ana[order((DE_ana$log2FoldChange), decreasing = TRUE), ]
DE_ana <- DE_ana[DE_ana$log2FoldChange < 0, ]
  
DE_spi <- res_spi_vs_other[complete.cases(res_spi_vs_other), ]
DE_spi <- DE_spi[DE_spi$padj < 0.01,]
DE_spi <- DE_spi[order((DE_spi$log2FoldChange), decreasing = TRUE), ]
DE_spi <- DE_spi[DE_spi$log2FoldChange > 0, ]
    
DE <- c(rownames(DE_spi),rownames(DE_salt), rownames(DE_ana))
  
heatmap_cells <- c(rownames(single[single$condition == "spi",]), rownames(single[single$condition == "salt",]), rownames(single[single$condition == "ana",]))  
norm <- nSingle[DE,heatmap_cells]
scaled <- t(apply(norm, 1, scale))
    
group <- factor(c(rep("spi", 19), rep("salt", 23), rep("ana", 27)))
ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("spi" = "#619cff", "salt" = "#00BA38", "ana" = "#F8766D")))
pdf("./output/heatmap_single_DEseq.pdf", width = 14, height = 7)
p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), top_annotation = ha, row_names_gp = gpar(fontsize = 2))
print(p)
dev.off()



write.csv(DE_salt, "./output/DE_salt_single.csv")
write.csv(DE_ana, "./output/DE_ana_single.csv")
write.csv(DE_spi, "./output/DE_spi_single.csv") 

################################################################################
################## Comparision to Jay Hinton paper --- Single ##################
################################################################################
DE_salt <- rownames(read.csv("./output/DE_salt_single.csv", row.names = 1))
DE_table_salt <- nSingle[DE_salt, rownames(single[single$condition == "salt",])]
DE_salt <- gsub("1344_", "", DE_salt)
rownames(DE_table_salt) <- DE_salt
DE_table_salt <- data.frame(rowMeans(DE_table_salt))
    
DE_ana <- rownames(read.csv("./output/DE_ana_single.csv", row.names = 1))
DE_table_ana <- nSingle[DE_ana, rownames(single[single$condition == "ana",])]
DE_ana <- gsub("1344_", "", DE_ana)
rownames(DE_table_ana) <- DE_ana
DE_table_ana <- data.frame(rowMeans(DE_table_ana))

tpm <- read_excel("./Kroger/mmc3.xlsx")
colnames(tpm) <- tpm[2,]
tpm <- tpm[-c(1,2), ]
tpm <- as.data.frame(tpm)
rownames(tpm) <- tpm[,1]
    
DE_salt <- intersect(DE_salt, rownames(tpm))
DE_table_salt <- DE_table_salt[DE_salt,]
tpm_salt <- tpm[DE_salt,c(19,22)]
tpm_salt$ratio <- as.numeric(tpm_salt$`NaCl shock`) / as.numeric(tpm_salt$`Anaerobic shock`)
tpm_salt$Group<-ifelse(tpm_salt$ratio>1,"A:Up in Kroger","B:Down in Kroger")
tpm_salt$Identity<-ifelse(tpm_salt$ratio>0,"A","B")
tpm_salt$gene <- rownames(tpm_salt)

## Fig2c_single  
ggplot(tpm_salt,aes(x=gene,y=ratio,fill=Group)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  scale_y_log10() +
  ggtitle("Up genes in Salt condition in our data (single)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  ylab("tpm salt / tpm ana (log scale)") +
  scale_fill_manual(values=c("#E69F00","#999999"))

ggsave("./output/Fig3c_single_salt_histogram.png")


ggplot(tpm_salt,aes(x=Identity,y=ratio)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_log10() +
  geom_hline(yintercept=1, linetype="dashed") 
  
ggsave("./output/Fig3c_single_salt_boxplot.png", width = 2)

   
    
DE_ana <- intersect(DE_ana, rownames(tpm))
DE_table_ana <- DE_table_ana[DE_ana,]
tpm_ana <- tpm[DE_ana,c(19,22)]
tpm_ana$ratio <- as.numeric(tpm_ana$`Anaerobic shock`) / as.numeric(tpm_ana$`NaCl shock`)
tpm_ana$Group<-ifelse(tpm_ana$ratio>1,"A: Up in Kroger","B: Down in Kroger")
tpm_ana$Identity<-ifelse(tpm_ana$ratio>0,"A","B")
tpm_ana$gene <- rownames(tpm_ana)
    

ggplot(tpm_ana,aes(x=gene,y=ratio,fill=Group)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  scale_y_log10() +
  ggtitle("Up genes in ana condition in our data (single)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  ylab("tpm ana / tpm salt (log scale)") +
  scale_fill_manual(values=c("#E69F00","#999999"))

ggsave("./output/Fig3c_single_salt_histogram.png")

ggplot(tpm_ana,aes(x=Identity,y=ratio)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_log10() +
  geom_hline(yintercept=1, linetype="dashed") 

ggsave("./output/Fig3c_single_salt_boxplot.png", width = 2)

#################################################################################
################## using DEseq2 - one vs other - ten cells ######################
#################################################################################
ten$SALT <- factor(as.character(gsub("spi", "other",gsub("ana", "other", ten$condition))))
ten$ANA <- factor(as.character(gsub("spi", "other",gsub("salt", "other", ten$condition))))
ten$SPI <- factor(as.character(gsub("salt", "other",gsub("ana", "other", ten$condition))))
    
#SALT
cts <- cbind(ten_table[,rownames(ten[ten$condition == "salt",])], ten_table[,rownames(ten[ten$condition == "ana",])], ten_table[,rownames(ten[ten$condition == "spi",])])
cts <- cts[rowSums(cts) > 0, ]
cts <- cts[rowSums(cts > 0) > 5, ]
  
col_data <- ten[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design= ~SALT)
dds <- DESeq(dds)
resultsNames(dds)
res_salt_vs_other <- results(dds, name = "SALT_salt_vs_other")
res_salt_vs_other <- data.frame(res_salt_vs_other)
res_salt_vs_other$gene_name <- rowData[rownames(res_salt_vs_other), ][,1] 
write.csv(res_salt_vs_other, file = "./output/res_salt_vs_other_ten.csv")  
    
    
#ANA
cts <- cbind(ten_table[,rownames(ten[ten$condition == "ana",])], ten_table[,rownames(ten[ten$condition == "salt",])], ten_table[,rownames(ten[ten$condition == "spi",])])
cts <- cts[rowSums(cts) > 0, ]
cts <- cts[rowSums(cts > 0) > 5, ]

col_data <- ten[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design= ~ANA)
dds <- DESeq(dds)
resultsNames(dds)
res_ana_vs_other <- results(dds, name = "ANA_other_vs_ana")
res_ana_vs_other <- data.frame(res_ana_vs_other)
res_ana_vs_other$gene_name <- rowData[rownames(res_ana_vs_other), ][,1] 
write.csv(res_ana_vs_other, file = "./output/res_ana_vs_other_ten.csv")  
    
    
    
#SPI
cts <- cbind(ten_table[,rownames(ten[ten$condition == "spi",])], ten_table[,rownames(ten[ten$condition == "salt",])], ten_table[,rownames(ten[ten$condition == "ana",])])
cts <- cts[rowSums(cts) > 0, ]
cts <- cts[rowSums(cts > 0) > 5, ]
  
col_data <- ten[colnames(cts),]
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design= ~SPI)
dds <- DESeq(dds)
resultsNames(dds)
res_spi_vs_other <- results(dds, name = "SPI_spi_vs_other")
res_spi_vs_other <- data.frame(res_spi_vs_other)
res_spi_vs_other$gene_name <- rowData[rownames(res_spi_vs_other), ][,1] 
write.csv(res_spi_vs_other, file = "./output/res_spi_vs_other_ten.csv")
    
    
    
    
DE_salt <- res_salt_vs_other[complete.cases(res_salt_vs_other), ]
DE_salt <- DE_salt[DE_salt$padj < 0.01,]
DE_salt <- DE_salt[order((DE_salt$log2FoldChange), decreasing = TRUE), ]
DE_salt <- DE_salt[DE_salt$log2FoldChange > 0, ]
  
DE_ana <- res_ana_vs_other[complete.cases(res_ana_vs_other), ]
DE_ana <- DE_ana[DE_ana$padj < 0.01,]
DE_ana <- DE_ana[order((DE_ana$log2FoldChange), decreasing = TRUE), ]
DE_ana <- DE_ana[DE_ana$log2FoldChange < 0, ]
    
DE_spi <- res_spi_vs_other[complete.cases(res_spi_vs_other), ]
DE_spi <- DE_spi[DE_spi$padj < 0.01,]
DE_spi <- DE_spi[order((DE_spi$log2FoldChange), decreasing = TRUE), ]
DE_spi <- DE_spi[DE_spi$log2FoldChange > 0, ]

DE <- c(rownames(DE_spi),rownames(DE_salt), rownames(DE_ana))
  
    
    
heatmap_cells <- c(rownames(ten[ten$condition == "spi",]), rownames(ten[ten$condition == "salt",]), rownames(ten[ten$condition == "ana",]))  
norm <- nTen[DE,heatmap_cells]
scaled <- t(apply(norm, 1, scale))
  
group <- factor(c(rep("spi", 20), rep("salt", 19), rep("ana", 18)))
ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("spi" = "#619cff", "salt" = "#00BA38", "ana" = "#F8766D")))
pdf("./output/heatmap_ten_DEseq.pdf", width = 14, height = 7)
p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), top_annotation = ha, row_names_gp = gpar(fontsize = 1))
print(p)
dev.off()
  
  
    
write.csv(DE_salt, "./output/DE_salt_ten.csv")
write.csv(DE_ana, "./output/DE_ana_ten.csv")
write.csv(DE_spi, "./output/DE_spi_ten.csv")

################################################################################
################## Comparision to Jay Hinton paper --- Ten ##################
################################################################################
DE_salt <- rownames(read.csv("./output/DE_salt_ten.csv", row.names = 1))
DE_table_salt <- nTen[DE_salt, rownames(ten[ten$condition == "salt",])]
DE_salt <- gsub("1344_", "", DE_salt)
rownames(DE_table_salt) <- DE_salt
DE_table_salt <- data.frame(rowMeans(DE_table_salt))

DE_ana <- rownames(read.csv("./output/DE_ana_ten.csv", row.names = 1))
DE_table_ana <- nTen[DE_ana, rownames(ten[ten$condition == "ana",])]
DE_ana <- gsub("1344_", "", DE_ana)
rownames(DE_table_ana) <- DE_ana
DE_table_ana <- data.frame(rowMeans(DE_table_ana))

tpm <- read_excel("./Kroger/mmc3.xlsx")
colnames(tpm) <- tpm[2,]
tpm <- tpm[-c(1,2), ]
tpm <- as.data.frame(tpm)
rownames(tpm) <- tpm[,1]

DE_salt <- intersect(DE_salt, rownames(tpm))
DE_table_salt <- DE_table_salt[DE_salt,]
tpm_salt <- tpm[DE_salt,c(19,22)]
tpm_salt$ratio <- as.numeric(tpm_salt$`NaCl shock`) / as.numeric(tpm_salt$`Anaerobic shock`)
tpm_salt$Group<-ifelse(tpm_salt$ratio>1,"A:Up in Kroger","B:Down in Kroger")
tpm_salt$Identity<-ifelse(tpm_salt$ratio>=0,"A","B")
tpm_salt$gene <- rownames(tpm_salt)

ggplot(tpm_salt,aes(x=gene,y=ratio,fill=Group)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  scale_y_log10() +
  ggtitle("Up genes in Salt condition in our data (ten)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  ylab("tpm salt / tpm ana (log scale)") +
  scale_fill_manual(values=c("#E69F00","#999999"))

ggsave("./output/Fig3c_ten_salt_histogram.png")

ggplot(tpm_salt,aes(x=Identity,y=ratio)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_log10() +
  geom_hline(yintercept=1, linetype="dashed") 

ggsave("./output/Fig3c_ten_salt_boxplot.png", width = 2)


DE_ana <- intersect(DE_ana, rownames(tpm))
DE_table_ana <- DE_table_ana[DE_ana,]
tpm_ana <- tpm[DE_ana,c(19,22)]
tpm_ana$ratio <- as.numeric(tpm_ana$`Anaerobic shock`) / as.numeric(tpm_ana$`NaCl shock`)
tpm_ana$Group<-ifelse(tpm_ana$ratio>1,"A: Up in Kroger","B: Down in Kroger")
tpm_ana$Identity<-ifelse(tpm_ana$ratio>0,"A","B")
tpm_ana$gene <- rownames(tpm_ana)


ggplot(tpm_ana,aes(x=gene,y=ratio,fill=Group)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  scale_y_log10() +
  ggtitle("Up genes in ana condition in our data (ten)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  ylab("tpm ana / tpm salt (log scale)") +
  scale_fill_manual(values=c("#E69F00","#999999"))

ggsave("./output/Fig3c_ten_ana_histogram.png")

ggplot(tpm_ana,aes(x=Identity,y=ratio)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_log10() +
  geom_hline(yintercept=1, linetype="dashed") 

ggsave("./output/Fig3c_ten_ana_boxplot.png", width = 2)

