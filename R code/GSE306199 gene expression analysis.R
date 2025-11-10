# load required packages
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
------------------------------------------------------

  
# loading raw data
raw <- read.csv("gene expression analysis/GSE306199_raw_counts_24hr.csv.gz")

# count matrix
count_data <- raw[, -1]
rownames(count_data) <- raw[, 1]
write.csv(count_data, "gene expression analysis/data/count_matrix.csv")

# sample metadata
colData <- data.frame(
  row.names = colnames(count_data),
  infection = factor(rep(c("Uninfected", "Infected"), each = 10)),
  treatment = factor(rep(c("Vehicle", "Dillapiole"), times = 2, each = 5))
)
write.csv(colData, "gene expression analysis/data/sample_metadata.csv")

# testing if row names and column names(samples) match
all(colnames(count_data) %in% rownames(colData))
all(colnames(count_data) == rownames(colData))

# set reference levels essential for comparisons
colData$infection <- relevel(colData$infection, ref = 'Uninfected')
colData$treatment <- relevel(colData$treatment, ref = 'Vehicle')
colData$group <- factor(paste(colData$infection, colData$treatment, sep = '_'))

# creating DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data, 
                       colData = colData,
                       design = ~ group)

# selecting genes with more than 10 total counts 
keeps <- rowSums(counts(dds)) >= 10
dds <- dds[keeps, ]

# the actual differential gene expression analysis
dds <- DESeq(dds)

#library sizes quality control
library_sizes <- colSums(counts(dds, normalized = TRUE), na.rm = TRUE)
cat('library size summary:\n')
print(summary(library_sizes))
png('sample_QC_library_sizes.png', width = 1200, height = 800, res = 120)
par(mar = c(12, 5, 4, 2) + 0.1)

barplot(
  library_sizes,
  main = 'QC: normalized library sizes per sample',
  ylab = 'total normalized counts',
  las = 2,
  col = 'red',
  border = NA,
  cex.names = 0.7
)
abline(h = pretty(library_sizes), col = 'gray80', lty = 'dotted')
dev.off()
# designing a function to automate file-saving process for significant genes (0.05 < adjusted p-value and log2FoldChange > 1) 
save_results <- function(res, name) {
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj), ]
  write.csv(df,here('results', 'DE', paste0(name, '_all_genes.csv')))
  
  sig <- subset(df, padj < 0.05 & abs(log2FoldChange) > 1)
  write.csv(sig, here('results', 'DE', paste0(name, '_sig_genes.csv')))
  
  top_up <- head(sig[order(sig$log2FoldChange, decreasing = TRUE), ], 10)
  top_down <- head(sig[order(sig$log2FoldChange), ], 10)
  
  write.csv(top_up, here('results', 'DE', paste0(name, '_top10_up.csv')))
  write.csv(top_down, here('results', 'DE', paste0(name, '_top10_down.csv')))
  
  return(list(sig = sig, top = list(up = top_up, down = top_down)))
}

# a function that creates heatmaps for top DE genes
top10_heatmap <- function(top_list, vsd, ann, filename, title) {
  genes <- c(rownames(top_list$up), rownames(top_list$down))
  
  if (length(genes) >= 2) {
    pheatmap(
      assay(vsd)[genes, ],
      annotation_col = ann,
      show_rownames = TRUE,
      main = title,
      fontsize_row = 8,
      filename = filename,
      width = 10, 
      height = 8
    )
  }
}
  

# function for MA plots 
ma_plot <- function(res, title, filename) {
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj), ]
  
  ggplot(df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = c('grey', 'red')) +
    theme_minimal() +
    labs(title = title,
         x = 'Mean expression',
         y = 'Log2 fold change') +
    theme(legend.position = 'none')
  
  ggsave(filename, width = 8, height = 6)
}
------------------------------------------------------
  
  
# differential expression comparisons
vehicle_effect <- results(dds, contrast = c('group', 'Infected_Vehicle', 'Uninfected_Vehicle'))
dillapiole_effect <- results(dds, contrast = c('group', 'Infected_Dillapiole', 'Uninfected_Dillapiole'))
uninf_treat_effect <- results(dds, contrast = c('group', 'Uninfected_Dillapiole', 'Uninfected_Vehicle'))
inf_treat_effect <- results(dds, contrast = c('group', 'Infected_Dillapiole', 'Infected_Vehicle'))


------------------------------------------------------

# Gene Annotation
## a function to add gene annotation and description columns to the results
add_gene_names <- function(gene_list) {
  df <- as.data.frame(gene_list)
  df$ensembl_gene_id <- rownames(df)
  
  ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  
  gene_names <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
    filters = 'ensembl_gene_id',
    values = df$ensembl_gene_id,
    mart = ensembl
  )
  
  gene_names$description <- gsub('\\[Source.*', '', gene_names$description)
  gene_names$description <- trimws(gene_names$description)
  
  merged <- merge(df, gene_names, by = 'ensembl_gene_id', all.x = TRUE)
  merged <- merged[, c('ensembl_gene_id', 'external_gene_name', 'description', 
                       setdiff(names(merged), c('ensembl_gene_id', 'external_gene_name', 'description')))]
  
  return(merged)
}

vehicle_effect_saved <- save_results(vehicle_effect, 'infected_vs_uninfected_vehicle')
dillapiole_effect_saved    <- save_results(dillapiole_effect, 'infected_vs_uninfected_dillapiole')
uninf_treat_effect_saved   <- save_results(uninf_treat_effect, 'dillapiole_vs_vehicle_uninfected')
inf_treat_effect_saved     <- save_results(inf_treat_effect, 'dillapiole_vs_vehicle_infected')

# annotating top 10 up and down regulated genes
vehicle_up_annotated   <- add_gene_names(vehicle_effect_saved$top$up)
vehicle_down_annotated <- add_gene_names(vehicle_effect_saved$top$down)
dillapiole_up_annotated   <- add_gene_names(dillapiole_effect_saved$top$up)
dillapiole_down_annotated <- add_gene_names(dillapiole_effect_saved$top$down)


write.csv(vehicle_up_annotated, here('results', 'annotated', 'infection_vehicle_top10_up_annotated.csv'))
write.csv(vehicle_down_annotated, here('results', 'annotated', 'infection_vehicle_top10_down_annotated.csv'))
write.csv(dillapiole_up_annotated, here('results', 'annotated', 'infection_dillapiole_top10_up_annotated.csv'))
write.csv(dillapiole_down_annotated, here('results', 'annotated', 'infection_dillapiole_top10_down_annotated.csv'))
------------------------------------------------------

  
# applying the ma-plot function on each comparison
ma_plot(vehicle_effect, 'MA Plot: Infected vs Uninfected (Vehicle)', 'infected_vs_uninfected_vehicle_MA.png')
ma_plot(dillapiole_effect, 'MA Plot: Infected vs Uninfected (Dillapiole)', 'infected_vs_uninfected_dillapiole_MA.png')
ma_plot(uninf_treat_effect, 'MA Plot: Dillapiole vs Vehicle (Uninfected)', 'dillapiole_vs_vehicle_uninfected_MA.png')
ma_plot(inf_treat_effect, 'MA Plot: Dillapiole vs Vehicle (Infected)', 'dillapiole_vs_vehicle_infected_MA.png')

# drawing pca plots to see the overall structure of the data and visualize sample relationships
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c('infection', 'treatment'), returnData = TRUE)

ggplot(pcaData, aes(PC1, PC2, color = infection, shape = treatment)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = 'PCA of All Samples')
ggsave('PCA_plot.png', width = 8, height = 6)
------------------------------------------------------

  
# applying top10_heatmap 
ann <- as.data.frame(colData(vsd)[, c('infection','treatment')])

top10_heatmap(vehicle_effect_saved$top, vsd, ann, 
                   'infected_vs_uninfected_vehicle_heatmap.png',
                   'Top 10 DE Genes - Inf and Uninf (Vehicle)')

top10_heatmap(dillapiole_effect_saved$top, vsd, ann, 
                   'infected_vs_uninfected_dillapiole_heatmap.png',
                   'Top 10 DE Genes - Inf and Uninf (Dillapiole)')

top10_heatmap(uninf_treat_effect_saved$top, vsd, ann, 
                   'dillapiole_vs_vehicle_uninfected_heatmap.png',
                   'Top 10 DE Genes - Treatment Effect (Uninfected)')

top10_heatmap(inf_treat_effect_saved$top, vsd, ann, 
                   'dillapiole_vs_vehicle_infected_heatmap.png',
                   'Top 10 DE Genes - Treatment Effect (Infected)')


# creating summary heatmaps combining multiple comparisons
infection_genes <- unique(c(
  rownames(vehicle_effect_saved$top$up), rownames(vehicle_effect_saved$top$down),
  rownames(dillapiole_effect_saved$top$up), rownames(dillapiole_effect_saved$top$down)
))

treatment_genes <- unique(c(
  rownames(uninf_treat_effect_saved$top$up), rownames(uninf_treat_effect_saved$top$down),
  rownames(inf_treat_effect_saved$top$up), rownames(inf_treat_effect_saved$top$down)
))


# infection heatmap
if (length(infection_genes) >= 2) {
  pheatmap(
    assay(vsd)[infection_genes, ],
    annotation_col = ann,
    main = 'Top Infection Genes',
    show_rownames = TRUE,
    fontsize_row = 8,
    filename = 'infection_heatmap.png',
    width = 10,
    height = 8
  )
}

# treatment heatmap
if (length(treatment_genes) >= 2) {
  pheatmap(
    assay(vsd)[treatment_genes, ],
    annotation_col = ann, 
    main = 'Top Treatment Genes',
    show_rownames = TRUE,
    fontsize_row = 8,
    filename = 'treatment_heatmap.png',
    width = 10,
    height = 8
  )
}
------------------------------------------------------
  
  
# volcano plots
EnhancedVolcano(as.data.frame(vehicle_effect),
                lab = rownames(vehicle_effect),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Infected vs Uninfected (Vehicle)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('infected_vs_uninfected_vehicle_volcano.png', width = 9, height = 7)

EnhancedVolcano(as.data.frame(dillapiole_effect),
                lab = rownames(dillapiole_effect),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Infected vs Uninfected (Dillapiole)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('infected_vs_uninfected_dillapiole_volcano.png', width = 9, height = 7)

EnhancedVolcano(as.data.frame(uninf_treat_effect),
                lab = rownames(uninf_treat_effect),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Dillapiole vs Vehicle (Uninfected)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('dillapiole_vs_vehicle_uninfected_volcano.png', width = 9, height = 7)

EnhancedVolcano(as.data.frame(inf_treat_effect),
                lab = rownames(inf_treat_effect),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Dillapiole vs Vehicle (Infected)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('dillapiole_vs_vehicle_infected_volcano.png', width = 9, height = 7)

