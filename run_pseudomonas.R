  library(ggplot2)
  library(gridExtra)
  library(DESeq2)
  library(scde)
  library(ComplexHeatmap)
  library(viridis)
  library(circlize)
  
  
  
  setwd("/address/to/working/directory/")
  
  ## importing row data
  anno <- read.table("./index/Pseudomonas_aeruginosa_pao1.ASM676v1.45.gtf", sep = "\t")
  anno <- anno[anno$V3 == "gene", ]
  
  biotype <- gsub(";", "", gsub(".*; gene_biotype ", "", anno$V9))
  gene_name <- gsub("; .*", "", gsub("gene_id.*", "",gsub(".*; gene_name ", "", anno$V9)))
  gene_id <- gsub("; .*", "",gsub("gene_id ", "", anno$V9))
  
  rowData <- data.frame(gene_name, biotype)
  rownames(rowData) <- gene_id
  
  ## Importing the count table
  temp <- list.files("./unique_count", pattern = "unique.txt")
  col <- gsub("Aligned.out.unique.txt", "", temp)
  filenames <- list.files(path = "./unique_count", full.names = TRUE)
  datalist <- lapply(filenames, function(x){read.table(file = x, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)[,6]})
  unique_table <- sapply(datalist, cbind)
  colnames(unique_table) <- col
  rownames(unique_table) <- rownames(rowData)
  
  
  ## colData 
  colData <- data.frame(library = col, n.cells = gsub(".*_", "", col)) 
  rownames(colData) <- col
  
  
  ## Adding the number of detected genes and library size
  colData$n.gene <- colSums(unique_table > 5)
  colData$lib.size <- colSums(unique_table)
  
  
  ## sFig6a
  p1 <- ggplot(data = colData, aes(x = n.gene, y = lib.size, color = n.cells)) +
    geom_point(size = 2) +
    geom_point(shape = 1,size = 2,colour = "black") +
    theme_classic() +
    theme(legend.position = "none") +
    ylab("Library size") +
    xlab("Number of detected genes")
  
  p2 <- ggplot(data = colData, aes(x = n.cells, y = n.gene, color = n.cells)) +
    geom_violin(alpha = 0) +
    geom_jitter(aes(fill=n.cells), colour="black",pch=21, size=2) +
    theme_classic() +
    theme(legend.position = "none") +
    ylab("Number of detected genes") +
    xlab("")
  
  
  sFig6a <- grid.arrange(p1, p2, ncol=2)
  ggsave("sFig6a.pdf", sFig6a, width = 7)
  
  ## sFig6b
  single <- colData[colData$n.cells == "sc",]
  ten <- colData[colData$n.cells == "10",]
  
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
    antisense = non_unique_table_single[grep("antisense", biotype), ],
    ribozyme = non_unique_table_single[grep("ribozyme", biotype), ],
    tmRNA = non_unique_table_single[grep("tmRNA", biotype), ]
  )
  
  percentage_single <- data.frame(
    rRNA = (colSums(non_unique_table_single[grep("rRNA", biotype), ]) / rowSums(ratio_single)) * 100,
    tRNA = (colSums(non_unique_table_single[grep("tRNA", biotype), ]) / rowSums(ratio_single)) * 100,
    ncRNA = (colSums(non_unique_table_single[grep("ncRNA", biotype), ]) / rowSums(ratio_single)) * 100,
    sRNA = (colSums(non_unique_table_single[grep("sRNA", biotype), ])  / rowSums(ratio_single)) * 100,
    protein_coding = (colSums(non_unique_table_single[grep("protein_coding", biotype), ])  / rowSums(ratio_single)) * 100,
    antisense = (non_unique_table_single[grep("antisense", biotype), ]  / rowSums(ratio_single)) * 100,
    ribozyme = (non_unique_table_single[grep("ribozyme", biotype), ]  / rowSums(ratio_single)) * 100,
    tmRNA = (non_unique_table_single[grep("tmRNA", biotype), ]  / rowSums(ratio_single)) * 100
  )
  
  write.csv(percentage_single, file = "percentage_single.csv")
  
  ## ten cells
  
  
  ratio_ten <- data.frame(
    rRNA = colSums(non_unique_table_ten[grep("rRNA", biotype), ]),
    tRNA = colSums(non_unique_table_ten[grep("tRNA", biotype), ]),
    ncRNA = colSums(non_unique_table_ten[grep("ncRNA", biotype), ]),
    sRNA = colSums(non_unique_table_ten[grep("sRNA", biotype), ]),
    protein_coding = colSums(non_unique_table_ten[grep("protein_coding", biotype), ]),
    antisense = non_unique_table_ten[grep("antisense", biotype), ],
    ribozyme = non_unique_table_ten[grep("ribozyme", biotype), ],
    tmRNA = non_unique_table_ten[grep("tmRNA", biotype), ]
  )
  
  percentage_ten <- data.frame(
    rRNA = (colSums(non_unique_table_ten[grep("rRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
    tRNA = (colSums(non_unique_table_ten[grep("tRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
    ncRNA = (colSums(non_unique_table_ten[grep("ncRNA", biotype), ]) / rowSums(ratio_ten)) * 100,
    sRNA = (colSums(non_unique_table_ten[grep("sRNA", biotype), ])  / rowSums(ratio_ten)) * 100,
    protein_coding = (colSums(non_unique_table_ten[grep("protein_coding", biotype), ])  / rowSums(ratio_ten)) * 100,
    antisense = (non_unique_table_ten[grep("antisense", biotype), ]  / rowSums(ratio_ten)) * 100,
    ribozyme = (non_unique_table_ten[grep("ribozyme", biotype), ]  / rowSums(ratio_ten)) * 100,
    tmRNA = (non_unique_table_ten[grep("tmRNA", biotype), ]  / rowSums(ratio_ten)) * 100
  )
  
  write.csv(percentage_ten, file = "percentage_ten.csv")
