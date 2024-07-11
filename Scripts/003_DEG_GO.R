
# DEG for comparing groups 
# GO with the results
# sel_sce is a SCE

deg_go_comparison <- function(sel_sce = qgo, group1 = "Up", group2 = "Down", grby = "TcrF_MM_BIC_Effect_050", 
                           folder_main = results_folder_full, folder = "DEG_GO")
{
  t1 <- Sys.time()
  print(t1)
  
  newdir <- paste0(folder_main, "/", folder)
  dir.create(newdir)
  setwd(newdir)
  
  df1 <- data.frame(table(colData(sel_sce)[,grby]))
  
  sel_sce$comparison_group <- ifelse(colData(sel_sce)[,grby] == group1, group1, 
                                     ifelse(colData(sel_sce)[,grby] == group2, group2, NA))
  sce_filtered <- sel_sce[, !is.na(sel_sce$comparison_group)]
  
  df2 <- data.frame(table(colData(sce_filtered)[,"comparison_group"]))
  
  # mean(apply(assay(qgo),1,mean))
  
  # dds <- DESeqDataSet(sel_sce)
  
  dds <- DESeqDataSet(sce_filtered, design = ~ comparison_group)
  dds <- DESeq(dds)
  
  t2<- Sys.time()
  print("Deseq")
  print(t1)
  print(t2)
  print(t2-t1)
  
  # resultsNames(dds)
  res <- results(dds, contrast = c("comparison_group", group1, group2))
  
  res <- results(dds)
  resData <- as.data.frame(res)
  
  resData$significativo <- resData$pvalue < 0.05  # Ajusta el umbral de p-valor según sea necesario
  resData$gen <- rownames(resData)
  
  head(resData)
  
  resData <- merge(resData, genedesc2, by.x = "gen", by.y = "GeneName", all.x = T)
  resData <- resData[order(resData$pvalue), ]
  
  # Hacemos el volcano plot
  gg1 <- ggplot(resData, aes(x = log2FoldChange, y = -log10(pvalue), color = significativo)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(x = "Log2 Fold Change", y = "-Log10 P-Valor") +
    geom_text_repel(
      aes(label = ifelse(resData$significativo & abs(resData$log2FoldChange) > 1, as.character(resData$gen), '')),
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = 'grey50'
    ) +
    scale_color_manual(values = c("grey", "red"))
  ggsave(paste0("Image_001_plot_clusters_",sample_name,"_", grby, "_", group1, "_", group2, "_ALL.jpg"), gg1, width = 15, height = 15)
  
  
  removegens <- c(grep("Ig", resData$gen),  grep("mt-", resData$gen), grep("Rps|Rpl", resData$gen))
  
  head(resData$gen, 20)
  resData2 <- resData[-removegens,]
  
  # Hacemos el volcano plot
  gg2 <- ggplot(resData2, aes(x = log2FoldChange, y = -log10(pvalue), color = significativo)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(x = "Log2 Fold Change", y = "-Log10 P-Valor") +
    geom_text_repel(
      aes(label = ifelse(resData2$significativo & abs(resData2$log2FoldChange) > 1, as.character(resData2$gen), '')),
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = 'grey50'
    ) +
    scale_color_manual(values = c("grey", "red"))
  ggsave(paste0("Image_001_plot_clusters_",sample_name,"_", grby, "_", group1, "_", group2, "_Filtered.jpg"), gg2, width = 15, height = 15)
  
  # Boxplot
  gene_data <- data.frame(expression = assay(sce_filtered)[resData2$gen[1], ], 
                          group = sce_filtered$comparison_group)
  
  gg3 <- ggplot(gene_data, aes(x = group, y = expression)) +
    geom_boxplot() +
    xlab("Grupo") +
    ylab(paste0("Expresión de ",resData2$gen[1])) +
    ggtitle(paste0("Comparación de la Expresión de ",resData2$gen[1] ," entre Grupos"))
  ggsave(paste0("Image_001_plot_clusters_",sample_name,"_", grby, "_", group1, "_", group2, "_boxplot_",resData2$gen[1],".jpg"), gg3, width = 15, height = 15)
  
  
  openxlsx::write.xlsx(resData, file = paste0("20240711_", grby, "_", group1, "_", group2, "_DEG_all.xlsx"))
  openxlsx::write.xlsx(resData2, file = paste0("20240711_", grby, "_", group1, "_", group2, "_DEG_filtrated.xlsx"))
  openxlsx::write.xlsx(df1, file = paste0("20240711_", grby, "_", group1, "_", group2, "_counts_all.xlsx"))
  openxlsx::write.xlsx(df2, file = paste0("20240711_", grby, "_", group1, "_", group2, "_counts_selected.xlsx"))
  
  # Pathways
  
  library(clusterProfiler)
  library(org.Mm.eg.db)
  
  # Filtrar genes significativamente diferentes
  resData <- resData[!is.na(resData$pvalue),]
  sig_genes <- resData$gen[resData$pvalue < 0.05]
  
  # Convertir nombres de genes a ENTREZ IDs
  entrez_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Análisis de enriquecimiento GO
  go_results_bp <- enrichGO(gene = entrez_ids$ENTREZID, 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "ENTREZID",
                            ont = "BP", # BP: Biological Process, MF: Molecular Function, CC: Cellular Component
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  
  # Visualizar los resultados
  b_BP <- barplot(go_results_bp, showCategory = 10)
  d_BP <- dotplot(go_results_bp, showCategory = 10)
  ggsave(paste0("Image_002_BP_barplot_",sample_name,"_", grby, "_", group1, "_", group2,".jpg"), b_BP, width = 15, height = 15)
  ggsave(paste0("Image_003_BP_dotplot_",sample_name,"_", grby, "_", group1, "_", group2,".jpg"), d_BP, width = 15, height = 15)
  
  go_results_MF <- enrichGO(gene = entrez_ids$ENTREZID, 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "ENTREZID",
                            ont = "MF", # BP: Biological Process, MF: Molecular Function, CC: Cellular Component
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  
  b_MF <- barplot(go_results_MF, showCategory = 10)
  d_MF <- dotplot(go_results_MF, showCategory = 10)
  ggsave(paste0("Image_002_MF_barplot_",sample_name,"_", grby, "_", group1, "_", group2,".jpg"), b_MF, width = 15, height = 15)
  ggsave(paste0("Image_003_MF_dotplot_",sample_name,"_", grby, "_", group1, "_", group2,".jpg"), d_MF, width = 15, height = 15)
  
  go_results_CC <- enrichGO(gene = entrez_ids$ENTREZID, 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "ENTREZID",
                            ont = "CC", # BP: Biological Process, MF: Molecular Function, CC: Cellular Component
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  
  b_CC <- barplot(go_results_CC, showCategory = 10)
  d_CC <- dotplot(go_results_CC, showCategory = 10)
  ggsave(paste0("Image_002_CC_barplot_",sample_name,"_", grby, "_", group1, "_", group2,".jpg"), b_CC, width = 15, height = 15)
  ggsave(paste0("Image_003_CC_dotplot_",sample_name,"_", grby, "_", group1, "_", group2,".jpg"), d_CC, width = 15, height = 15)
  
  # Export GO BP
  df_bp <- data.frame(go_results_bp)
  # Convertir los ENTREZ IDs de los resultados GO a símbolos de genes
  gene_symbols <- bitr(geneID = unlist(strsplit(df_bp$geneID, "/")), 
                       fromType = "ENTREZID", 
                       toType = "SYMBOL", 
                       OrgDb = org.Mm.eg.db)
  
  # Crear un vector para mapear los ENTREZ IDs a los símbolos de genes
  entrez_to_symbol <- setNames(gene_symbols$SYMBOL, gene_symbols$ENTREZID)
  
  # Función para reemplazar los ENTREZ IDs por símbolos de genes
  replace_entrez_with_symbol <- function(entrez_ids) {
    symbols <- entrez_to_symbol[unlist(strsplit(entrez_ids, "/"))]
    paste(symbols, collapse = "/")
  }
  
  # Aplicar la función a la columna geneID de los resultados GO
  df_bp$geneSymbol <- sapply(df_bp$geneID, replace_entrez_with_symbol)
  
  # Visualizar los resultados actualizados
  head(df_bp)
  openxlsx::write.xlsx(df_bp, file = paste0("20240711_", grby, "_", group1, "_", group2, "_GO_BP.xlsx"))
  
  # Export MF
  df_MF <- data.frame(go_results_MF)
  # Convertir los ENTREZ IDs de los resultados GO a símbolos de genes
  gene_symbols <- bitr(geneID = unlist(strsplit(df_MF$geneID, "/")), 
                       fromType = "ENTREZID", 
                       toType = "SYMBOL", 
                       OrgDb = org.Mm.eg.db)
  
  # Crear un vector para mapear los ENTREZ IDs a los símbolos de genes
  entrez_to_symbol <- setNames(gene_symbols$SYMBOL, gene_symbols$ENTREZID)
  
  # Aplicar la función a la columna geneID de los resultados GO
  df_MF$geneSymbol <- sapply(df_MF$geneID, replace_entrez_with_symbol)
  
  # Visualizar los resultados actualizados
  head(df_MF)
  openxlsx::write.xlsx(df_MF, file = paste0("20240711_", grby, "_", group1, "_", group2, "_GO_MF.xlsx"))
  
  # Export CC
  df_CC <- data.frame(go_results_CC)
  # Convertir los ENTREZ IDs de los resultados GO a símbolos de genes
  gene_symbols <- bitr(geneID = unlist(strsplit(df_CC$geneID, "/")), 
                       fromType = "ENTREZID", 
                       toType = "SYMBOL", 
                       OrgDb = org.Mm.eg.db)
  
  # Crear un vector para mapear los ENTREZ IDs a los símbolos de genes
  entrez_to_symbol <- setNames(gene_symbols$SYMBOL, gene_symbols$ENTREZID)
  
  # Aplicar la función a la columna geneID de los resultados GO
  df_CC$geneSymbol <- sapply(df_CC$geneID, replace_entrez_with_symbol)
  
  # Visualizar los resultados actualizados
  head(df_CC)
  openxlsx::write.xlsx(df_CC, file = paste0("20240711_", grby, "_", group1, "_", group2, "_GO_CC.xlsx"))
  
  t3<- Sys.time()
  print("End")
  print(t1)
  print(t2)
  print(t3)
  print(t3-t1)
  
}
