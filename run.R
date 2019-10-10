library(ggplot2)
library(gridExtra)
library(DESeq2)
library(scde)
library(ComplexHeatmap)
library(viridis)
library(circlize)


## importing col Data and removing one single outlier
colData <- read.csv("/home/ehsan/Desktop/salmonella/lib_changed_corr_fab.csv", row.names = 1)[-c(152,141), ]
## importing row data
anno <- read.table("/home/ehsan/Desktop/salmonella/index/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_sl1344.ASM21085v2.45.gtf", sep = "\t")
anno <- anno[anno$V3 == "gene", ]

biotype <- gsub(";", "", gsub(".*; gene_biotype ", "", anno$V9))
gene_name <- gsub("; .*", "", gsub("gene_id.*", "",gsub(".*; gene_name ", "", anno$V9)))
gene_id <- gsub("; .*", "",gsub("gene_id ", "", anno$V9))

rowData <- data.frame(gene_name, biotype)
rownames(rowData) <- gene_id

## Importing the count table
temp <- list.files("/home/ehsan/Desktop/salmonella/unique_count", pattern = "unique.txt")
col <- gsub(".*_", "", gsub("Aligned.sortedByCoord.out.unique.txt", "", temp))

filenames <- list.files(path = "/home/ehsan/Desktop/salmonella/unique_count", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,6]})
unique_table <- sapply(datalist, cbind)
colnames(unique_table) <- col
rownames(unique_table) <- rownames(read.table(file = "/home/ehsan/Desktop/salmonella/unique_count/L1800565_M1A1Aligned.sortedByCoord.out.unique.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE))


## removing the two outlier libraries
# colData <- colData[-c(152, 186),]

## limiting data to single and ten cells library
colData <- colData[colData$n.cells == "1" | colData$n.cells == "10",]
## similar order of libraries in colDAta and countTable 
unique_table <- unique_table[,rownames(colData)]

## Adding the number of detected genes and library size
colData$n.gene <- colSums(unique_table > 0)
colData$lib.size <- colSums(unique_table)

## Fig1c
p1 <- ggplot(data = colData[colData$n.cells == "1", ], aes(x = n.cells, y = n.gene, color = n.cells)) +
  geom_jitter() + 
  geom_violin(alpha = 0) +
  theme_classic() +
  ylab("Number of detected genes") +
  xlab("") +
  theme(legend.position = "none") +
  geom_hline(yintercept=50, linetype="dashed") +
  ggtitle("Single Cell")

p2 <- ggplot(data = colData[colData$n.cells == "10", ], aes(x = n.cells, y = n.gene, color = n.cells)) +
  geom_jitter() + 
  geom_violin(alpha = 0) +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  geom_hline(yintercept=50, linetype="dashed") +
  ggtitle("10 pooled cells")

Fig1c <- grid.arrange(p1, p2, ncol=2)
ggsave("Fig1c.svg", Fig1c)


## Fig1d
p1 <- ggplot(data = colData[colData$n.cells == "1", ], aes(x = condition, y=n.gene, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  geom_hline(yintercept=50, linetype="dashed") +
  ggtitle("Single Cell")


p2 <- ggplot(data = colData[colData$n.cells == "10", ], aes(x = condition, y=n.gene, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  geom_hline(yintercept=50, linetype="dashed") +
  ggtitle("10 pooled cells")

Fig1d <- grid.arrange(p1, p2, ncol=1)
ggsave("Fig1d.svg", Fig1d)


## sFig2b
p1 <- ggplot(data = colData[colData$n.cells == "1", ], aes(x = lib.size, y=n.gene, color = condition)) +
  geom_jitter() +
  theme_classic() +
  geom_hline(yintercept=50, linetype="dashed") +
  geom_vline(xintercept=5000, linetype="dashed") +
  theme(legend.position = "none") +
  ylab("Number of detected genes") +
  xlab("Library Size")
  
p2 <- ggplot(data = colData[colData$n.cells == "10", ], aes(x = lib.size, y=n.gene, color = condition)) +
  geom_jitter() +
  theme_classic() +
  geom_hline(yintercept=50, linetype="dashed") +
  geom_vline(xintercept=20000, linetype="dashed") +
  theme(legend.position = "none") +
  ylab("Number of detected genes") +
  xlab("Library Size")

sFig2b <- grid.arrange(p1, p2, ncol=2)
ggsave("sFig2b.svg", sFig2b, width = 12)


## sFig2a
p1 <- ggplot(data = colData[colData$n.cells == "1", ], aes(x = condition, y=lib.size, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  geom_hline(yintercept=5000, linetype="dashed") +
  ggtitle("Single Cell")


p2 <- ggplot(data = colData[colData$n.cells == "10", ], aes(x = condition, y=lib.size, color = condition)) +
  geom_jitter() +
  geom_violin(alpha = 0) + 
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none") + 
  geom_hline(yintercept=20000, linetype="dashed") +
  ggtitle("10 pooled cells")

sFig2a <- grid.arrange(p1, p2, ncol=2)
ggsave("sFig2a.svg", sFig2a, width = 14)



### Filtering low quality cells
single <- colData[colData$n.cells == "1" & colData$n.gene > 50 & colData$lib.size > 5000, ]
###### removing the single spi library
single <- single[-1,]


ten <- colData[colData$n.cells == "10" & colData$n.gene > 50 & colData$lib.size > 20000, ]

single_table <- unique_table[,rownames(single)]
ten_table <- unique_table[,rownames(ten)]


### Normalizing data
sfSingle <- estimateSizeFactorsForMatrix(single_table)
sfTen <- estimateSizeFactorsForMatrix(ten_table)

nSingle <- t(t(single_table) / sfSingle)
nTen <- t(t(ten_table) / sfTen)

### Selecting the top variable genes and log transformation
topVarGenesSingle <- head(order(rowVars(nSingle), decreasing = TRUE),500)
nSingle_fin <- log(nSingle[topVarGenesSingle,] + 1)

topVarGenesTen <- head(order(rowVars(nTen), decreasing = TRUE),500)
nTen_fin <- log(nTen[topVarGenesTen,] + 1)


### Performing PCA (ten cells)
res.pca <- prcomp(t(nTen_fin), scale = FALSE)
pca_ten <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = ten$condition), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") 
  
  

### Performing PCA (single cell)
res.pca <- prcomp(t(nSingle_fin), scale = FALSE)
pca_single <- ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], condition = single$condition), aes(x = PC1, y = PC2, color = condition)) +
  geom_jitter(size = 2) +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") 

Fig2a <- grid.arrange(pca_single, pca_ten, ncol=1)
ggsave("Fig2a.svg", Fig2a, height = 12)


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
ggsave("sFig4a.svg", sFig4a, width = 20)


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

sFig4b <- grid.arrange(plot_salt, plot_ana, ncol = 3)
ggsave("sFig4b.svg", sFig4b, width = 20)


## Fig1b
### Here the non_unique counts are used to calculate the percentage of rRNA, tRNA, etc.
filenames <- list.files(path = "/home/ehsan/Desktop/salmonella/non_unique_count", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,6]})
non_unique_table <- sapply(datalist, cbind)
colnames(non_unique_table) <- col
rownames(non_unique_table) <- rownames(rowData)


non_unique_table_single <- non_unique_table[,rownames(single)]
non_unique_table_ten <- non_unique_table[,rownames(ten)]

## calculating intergenic regions for all libraries
filenames <- list.files(path = "/home/ehsan/Desktop/salmonella/non_unique_summary", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,1]})
non_unique_summary <- sapply(datalist, cbind)
colnames(non_unique_summary) <- col
IGR <- t(non_unique_summary[10,])


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
  antitoxin = colSums(non_unique_table_single[grep("antitoxin", biotype), ]),
  IGR = IGR[,rownames(single)]
)

write.csv(ratio_single, file = "ratio_single.csv")


percentage_single <- data.frame(
  rRNA = (colSums(non_unique_table_single[grep("rRNA", biotype), ]) / rowSums(ratio_single)) * 100,
  tRNA = (colSums(non_unique_table_single[grep("tRNA", biotype), ]) / rowSums(ratio_single)) * 100,
  ncRNA = (colSums(non_unique_table_single[grep("ncRNA", biotype), ]) / rowSums(ratio_single)) * 100,
  sRNA = (colSums(non_unique_table_single[grep("sRNA", biotype), ])  / rowSums(ratio_single)) * 100,
  protein_coding = (colSums(non_unique_table_single[grep("protein_coding", biotype), ])  / rowSums(ratio_single)) * 100,
  pseudogene = (colSums(non_unique_table_single[grep("pseudogene", biotype), ])  / rowSums(ratio_single)) * 100,
  antisense = (colSums(non_unique_table_single[grep("antisense", biotype), ])  / rowSums(ratio_single)) * 100,
  ribozyme = (non_unique_table_single[grep("ribozyme", biotype), ]  / rowSums(ratio_single)) * 100,
  antitoxin = (colSums(non_unique_table_single[grep("antitoxin", biotype), ]) / rowSums(ratio_single)) * 100,
  IGR = (IGR[,rownames(single)]  / rowSums(ratio_single)) * 100
)

write.csv(percentage_single, file = "percentage_single.csv")
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
  antitoxin = colSums(non_unique_table_ten[grep("antitoxin", biotype), ]),
  IGR = IGR[,rownames(ten)]
)

write.csv(ratio_ten, file = "ratio_ten.csv")



percentage_ten <- data.frame(
  rRNA = (colSums(non_unique_table_ten[grep("rRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
  tRNA = (colSums(non_unique_table_ten[grep("tRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
  ncRNA = (colSums(non_unique_table_ten[grep("ncRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
  sRNA = (colSums(non_unique_table_ten[grep("sRNA", biotype), ])  / rowSums(ratio_ten)) * 100,
  protein_coding = (colSums(non_unique_table_ten[grep("protein_coding", biotype), ])  / rowSums(ratio_ten)) * 100,
  pseudogene = (colSums(non_unique_table_ten[grep("pseudogene", biotype), ])  / rowSums(ratio_ten)) * 100,
  antisense = (colSums(non_unique_table_ten[grep("antisense", biotype), ])  / rowSums(ratio_ten)) * 100,
  ribozyme = (non_unique_table_ten[grep("ribozyme", biotype), ]  / rowSums(ratio_ten)) * 100,
  antitoxin = (colSums(non_unique_table_ten[grep("antitoxin", biotype), ]) / rowSums(ratio_ten)) * 100,
  IGR = (IGR[,rownames(ten)]  / rowSums(ratio_ten)) * 100
)

write.csv(percentage_ten, file = "percentage_ten.csv")



## Running the differential expression analysis with scde
## Single cell SALT vs ANA
cd <- cbind(single_table[,rownames(single[single$condition == "salt",])], single_table[,rownames(single[single$condition == "ana",])])
group <- factor(c(rep("salt", 13), rep("ana", 20)))
names(group) <- colnames(cd)
o.ifm <- scde.error.models(counts = cd, groups = , n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
ediff <- ediff[order((ediff$Z), decreasing = TRUE), ]
ediff$gene_name <- rowData[rownames(ediff), ][,1]
write.csv(ediff, file = "SALTvsANA_single.csv")


## creating the heatmap for single cell
DE <- read.csv("SALTvsANA_single.csv", header = TRUE, row.names = 1)
DE <- DE[abs(DE$cZ) > 1.28, ]

#count <- unique_table[rownames(DE),c(rownames(colData[colData$condition == "spi" & colData$n.cells == 1,]))]
norm <- nSingle[rownames(DE),colnames(cd)]
#norm <- cbind(count, norm)
scaled <- t(apply(norm, 1, scale))

#group <- data.frame(group = c(rep("spi", 30), rep("salt", 13), rep("ana", 20)))
ha <- HeatmapAnnotation(df = data.frame(group), col = list(group = c("spi" = "#619cff","salt" = "#00BA38", "ana" = "#F8766D")))

pdf("heatmap_single.pdf", width = 14, height = 7)
p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), top_annotation = ha)
print(p)
dev.off()

## Ten cells SALT vs Other
cd <- cbind(ten_table[,rownames(ten[ten$condition == "salt",])], ten_table[,rownames(ten[ten$condition == "ana" | ten$condition == "spi",])])
group <- factor(c(rep("salt", 33), rep("other", 33)))
names(group) <- colnames(cd)
o.ifm <- scde.error.models(counts = cd, groups = , n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
ediff <- ediff[order((ediff$Z), decreasing = TRUE), ]
ediff$gene_name <- rowData[rownames(ediff), ][,1]
write.csv(ediff, file = "SALTvsOTHER_ten.csv")

## Ten cells ANA vs Other
cd <- cbind(ten_table[,rownames(ten[ten$condition == "ana",])], ten_table[,rownames(ten[ten$condition == "salt" | ten$condition == "spi",])])
group <- factor(c(rep("ana", 13), rep("other", 53)))
names(group) <- colnames(cd)
o.ifm <- scde.error.models(counts = cd, groups = , n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
ediff <- ediff[order((ediff$Z), decreasing = TRUE), ]
ediff$gene_name <- rowData[rownames(ediff), ][,1]
write.csv(ediff, file = "ANAvsOTHER_ten.csv")


## Ten cells SPI vs Other
cd <- cbind(ten_table[,rownames(ten[ten$condition == "spi",])], ten_table[,rownames(ten[ten$condition == "salt" | ten$condition == "ana",])])
group <- factor(c(rep("spi", 20), rep("other", 46)))
names(group) <- colnames(cd)
o.ifm <- scde.error.models(counts = cd, groups = , n.cores = 8, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  group, n.randomizations  =  100, n.cores  =  8, verbose  =  1)
ediff <- ediff[order((ediff$Z), decreasing = TRUE), ]
ediff$gene_name <- rowData[rownames(ediff), ][,1]
write.csv(ediff, file = "SPIvsOTHER_ten.csv")

## Creating heatmap for ten cells
DE_salt <- read.csv("SALTvsOTHER_ten.csv", header = TRUE, row.names = 1)
DE_salt <- DE_salt[DE_salt$cZ < -1.96,]
DE_ana <- read.csv("ANAvsOTHER_ten.csv", header = TRUE, row.names = 1)
DE_ana <- DE_ana[DE_ana$cZ > 1.96,]
DE_spi <- read.csv("SPIvsOTHER_ten.csv", header = TRUE, row.names = 1)
DE_spi <- DE_spi[DE_spi$cZ < -1.96,]

DE <- c(rownames(DE_spi), rownames(DE_salt), rownames(DE_ana))

norm <- nTen[DE, c(rownames(ten[ten$condition == "spi", ]), rownames(ten[ten$condition == "salt", ]), rownames(ten[ten$condition == "ana", ]))]
scaled <- t(apply(norm, 1, scale))
group <- data.frame(group = c(rep("spi", 20), rep("salt", 33), rep("ana", 13)))
ha <- HeatmapAnnotation(df = group, col = list(group = c("spi" = "#619cff" ,"salt" = "#00BA38", "ana" = "#F8766D")))

pdf("heatmap_ten.pdf", width = 14, height = 7)
p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), top_annotation = ha, row_names_gp = gpar(fontsize = 5))
print(p)
dev.off()

