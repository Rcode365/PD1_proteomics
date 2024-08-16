
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#################  1. install packages  ############### 
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

all_packages <- c("tidyverse", "stringr", "pheatmap", "limma", "devtools", "clusterProfiler", "genefilter", "matrixStats", "ConsensusClusterPlus", "corrplot", "PCAtools", "UniprotR", "GseaVis", "xlsx","openxlsx","survival","survminer","installr")
installed <- installed.packages()
no_install <- all_packages[!all_packages %in% installed]
BiocManager::install(no_install, dependencies = TRUE, ask = FALSE, force = TRUE)


#################  2. load packages   ###############   
for(i in all_packages) { 
  library(i, character.only = TRUE) 
}

#################  3. install proteoDA   ############### 
##In R GUI, perform the following commands:

install.Rtools(keep_install_file=TRUE)

devtools::install_github("ByrumLab/proteoDA", dependencies = TRUE, build_vignettes = TRUE)

library(proteoDA)

#################  4.	Read the proteomics results

ms <- read.delim("proteinGroups.txt", header = TRUE, sep = "\t", check.names = FALSE)
dim(ms)
#[1] 836 698

## or load the demo data, which includes protein identification results and metadata for sample grouping.
load("SupplementaryData.RData")

#################  5.	 Filtering of protein identifications.   
ms_LFQ <- ms[ms$`Only identified by site`!="+" & ms$Reverse!="+" & ms$`Potential contaminant`!="+" & ms$`Unique peptides`>=2, ]
ms_LFQ <- ms_LFQ %>% select('Majority protein IDs','Gene names','Fasta headers','Protein names',"LFQ intensity B039":"LFQ intensity N-9")

#################  6.	Change protein (row) names.   
uniprot_id <- str_split(ms_LFQ$'Fasta headers',";") %>% sapply("[",1)
Acc_number <- str_split(ms_LFQ$'Majority protein IDs',";") %>% sapply("[",1)
Gene_name <- strsplit(ms_LFQ$'Gene names',";") %>% sapply("[",1)
Protein_name <- ms_LFQ$`Protein names`
row.names(ms_LFQ) <- uniprot_id  
ms_LFQ[,c(1:4)] <- NULL

#################  7.	Generate a data frame of protein annotations.   
protein_anno <- data.frame(Acc_number, uniprot_id, Gene_name, Protein_name, row.names = uniprot_id)

## Several symbols are unavailable, substitute with UniProtI ID
NA_gene <- protein_anno$uniprot_id[which(is.na(protein_anno$Gene_name))] %>%
  strsplit("_") %>% sapply("[",1)
protein_anno$Gene_name[which(is.na(protein_anno$Gene_name))] <- NA_gene
xlsx::write.xlsx(protein_anno,"protein_anno.xlsx")


#################  8.	Change sample (column) names.   
sample_name <- str_split(colnames(ms_LFQ)," ") %>%  sapply("[", 3)
colnames(ms_LFQ) <- sample_name
ms_LFQ <- ms_LFQ[,order(names(ms_LFQ))]  ## adjust sample order
head(ms_LFQ)[1:5, 1:5]  ## list the final matrix



#################  9.	Generate a sample group annotation.   

gp <- read.csv("group.txt", header = T, sep = "\t")
row.names(gp) <- gp[,1]
gp[,3] <- NULL
gp <- gp[order(row.names(gp)), ]  ## adjust the sample order
head(gp)


#################  10.	Building a DAList 
DAList_nor <- DAList(data = ms_LFQ,
                     annotation = protein_anno,
                     metadata = gp,
                     design = NULL,
                     eBayes_fit = NULL,
                     results = NULL,
                     tags = NULL)

#################  11.	 Filter out unwanted proteins with too many missing values
DAList_nor_filtered <- DAList_nor %>% 
  zero_to_missing() %>%
  filter_proteins_by_proportion(min_prop = 0.1, grouping_column = "condition")


#################  12.	Make the normalization report
write_norm_report(DAList_nor_filtered,
                  grouping_column = "condition",
                  output_dir = "Normalization_report",
                  filename = NULL,  ## default to normalization_report.pdf
                  overwrite = TRUE,
                  suppress_zoom_legend = FALSE,
                  use_ggrastr = FALSE)

#################  13.	Normalize the protein intensity using the most suitable method
ms_LFQ_nor <- normalize_data(DAList_nor_filtered, norm_method = "vsn")




#################  14.	Make the quality control report
write_qc_report(
  ms_LFQ_nor,
  color_column = "condition",
  label_column = NULL,
  output_dir = "Normalization_report",
  filename = NULL,   ## default to "QC_Report.pdf"
  overwrite = TRUE,
  top_proteins = 500,
  standardize = TRUE,
  pca_axes = c(1, 2),  ## using the first two principle components for PCA analysis
  dist_metric = "euclidean",
  clust_method = "complete",
  show_all_proteins = FALSE
)



#################  15.	Convert sample group in metadata to a factor and define the desired levels.
ms_LFQ_nor$metadata$condition <- factor(ms_LFQ_nor$metadata$condition,
                             levels = c("Normal","Post_NR", "Post_R","Pre_NR", "Pre_R"))
ms_LFQ_nor$metadata$condition



#################  16.	Build a no-intercept model
no_intercept <- add_design(ms_LFQ_nor, design_formula = ~0 + condition)

#################  17.	Make contrasts
no_intercept <- add_contrasts(no_intercept, contrasts_vector = c("Pre_R_vs_Pre_NR = Pre_R - Pre_NR", "Pre_R_vs_Post_NR = Pre_R - Post_NR", "Pre_R_vs_Post_R = Pre_R - Post_R", "Post_R_vs_Post_NR = Post_R - Post_NR", "Post_R_vs_Pre_NR = Post_R - Pre_NR", "Pre_NR_vs_Post_NR = Pre_NR - Post_NR", "Pre_R_vs_Normal = Pre_R - Normal", "Pre_NR_vs_Normal = Pre_NR - Normal", "Post_NR_vs_Normal = Post_NR - Normal", "Post_R_vs_Normal = Post_R - Normal"))



#################  18.	Fit a limma model
fit <- fit_limma_model(no_intercept)


#################  19.	Extract the limma results
results <- extract_DA_results(fit,
                              pval_thresh = 0.05,
                              lfc_thresh = 0.25,
                              adj_method = "BH", ## Benjamini-Hochberg
                              extract_intercept = FALSE)



#################  20.	Export the limma results
write_limma_tables(results,
                   output_dir = "limma",
                   overwrite = TRUE,
                   contrasts_subdir = "per_contrast_results",
                   summary_csv = NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = TRUE)


#################  21.	Write the plots of limma results and generate a multi-sheet table 
write_limma_plots(results,
                  grouping_column = "condition",
                  table_columns = c("uniprot_id"),
                  title_column = NULL,
                  output_dir = "limma/limma_plots",
                  tmp_subdir = "tmp",
                  overwrite = TRUE,
                  height = 1000,
                  width = 1000)

tab <- openxlsx::createWorkbook()
for (contrast in no_intercept$design$contrast_vector) {
  contrast <- strsplit(contrast, " ") %>% sapply("[", 1)
  contrast_result <- results$results[[contrast]]
  contrast_result <- contrast_result[contrast_result$adj.P.Val < 0.05,]
  openxlsx::addWorksheet(tab, contrast)
  openxlsx::writeData(tab, contrast, contrast_result, rowNames = TRUE)
}
openxlsx::saveWorkbook(tab,"limma_allcontrast.xlsx", overwrite = TRUE)


#################  22.	Extract the expression matrix of the differential proteins between responders and non-responders.
attach(results$results$Pre_R_vs_Pre_NR)
R_vs_NR_diff <- results$results$Pre_R_vs_Pre_NR[adj.P.Val < 0.05 & abs(logFC) > 0.25, ]
detach(results$results$Pre_R_vs_Pre_NR)

R_vs_NR_diff_data <- ms_LFQ_nor$data[rownames(R_vs_NR_diff), ]

gp2 <- gp[order(gp$condition), ]
R_vs_NR_diff_data2 <- R_vs_NR_diff_data[, rownames(gp2)]


#################  23.	Hierarchical cluster analysis using R package pheatmap. 
pdf(file = "heatmap.pdf")
pheatmap(R_vs_NR_diff_data2,
          scale = "row",
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          color = colorRampPalette(colors = c("blue", "white", "red"))(100),
          fontsize = 6,
          cellwidth = 3,
          cellheight = 5,
          border_color = "NA",
          annotation_col = gp2[2],
          cutree_rows = 2)
dev.off()


#################  24.	Extract gene symbols, accessions, and names of each protein.
attach(results$results$Pre_R_vs_Pre_NR)
up_pro <- rownames(results$results$Pre_R_vs_Pre_NR[logFC > 0.25 & adj.P.Val < 0.05, ])
up_symbol <- protein_anno[protein_anno$uniprot_id %in% up_pro, ]$Gene_name
up_accession <- protein_anno[protein_anno$uniprot_id %in% up_pro, ]$Acc_number

down_pro <- rownames(results$results$Pre_R_vs_Pre_NR[logFC < 0.25 & adj.P.Val < 0.05, ])
down_symbol <- protein_anno[protein_anno$uniprot_id %in% down_pro, ]$Gene_name 
down_accession <- protein_anno[protein_anno$uniprot_id %in% down_pro, ]$Acc_number
detach(results$results$Pre_R_vs_Pre_NR)

length(up_symbol)  
length(down_symbol)  


#################  25.	Enrichment analysis using R package UniprotR.
Enrichment.BP(up_symbol, OS = "hsapiens", p_value = 0.05, top = 15)  
Enrichment.BP(down_symbol, OS = "hsapiens", p_value = 0.05, top = 15)
Enrichment.KEGG(up_symbol, OS = "hsapiens", p_value = 0.05, top = 15)
Enrichment.KEGG(down_symbol, OS = "hsapiens", p_value = 0.05, top = 15)
Enrichment.MF(up_symbol, OS = "hsapiens", p_value = 0.05, top = 15)
Enrichment.MF(down_symbol,  OS = "hsapiens", p_value = 0.05, top = 15)
Enrichment.REAC(up_symbol, OS = "hsapiens", p_value = 0.05, top = 15)
Enrichment.REAC(down_symbol, OS = "hsapiens", p_value = 0.05, top = 15)


#################  26.	Prepare an ordered list for the gene set enrichment analysis
attach(results$results$Pre_R_vs_Pre_NR)
Pre_R_vs_Pre_NR_order <- results$results$Pre_R_vs_Pre_NR[order(logFC, decreasing =T), ]
detach(results$results$Pre_R_vs_Pre_NR)

genelist <- rownames(Pre_R_vs_Pre_NR_order)
## extract table
protein_anno_R_NR <- protein_anno[genelist, ]

## produce a data frame
symbol_FC <- data.frame(Pre_R_vs_Pre_NR_order$logFC, 
                         row.names = protein_anno_R_NR$Gene_name)
colnames(symbol_FC) <- "logFC"
head(symbol_FC)

## convert the data frame to a name vector
symbol_FC_vector <- setNames(symbol_FC$logFC, rownames(symbol_FC))
head(symbol_FC_vector)


#################  27.	Gene set enrichment analysis against the KEGG pathway database
cp_kegg <- clusterProfiler::read.gmt("c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")

GSEA_kegg <- clusterProfiler::GSEA(symbol_FC_vector, 
                                    TERM2GENE = cp_kegg, 
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05)

pdf("GSEA_KEGG.pdf")
GSEA_kegg2 <- as.data.frame(GSEA_kegg)
GseaVis::gseaNb(object = GSEA_kegg, geneSetID = rownames(GSEA_kegg2)[1], addPval = TRUE)
dev.off()

#################  28.	Gene set enrichment analysis against the Wikipathways database
wikipathways <- clusterProfiler::read.gmt("c2.cp.wikipathways.v2023.2.Hs.symbols.gmt")

GSEA_wikipathways <- clusterProfiler::GSEA(symbol_FC_vector, 
                                            TERM2GENE = wikipathways, 
                                            pAdjustMethod = "BH", 
                                            pvalueCutoff = 0.05)

pdf("GSEA_wikipathways.pdf")
GSEA_wikipathways2 <- as.data.frame(GSEA_wikipathways)
GseaVis::gseaNb(object = GSEA_wikipathways,
        geneSetID = rownames(GSEA_wikipathways2)[1], addPval = TRUE)
dev.off()


#################  29.	Prepare the input matrix for PCA analysis
ms_LFQ_nor_data <- as.data.frame(ms_LFQ_nor$data)
ms_LFQ_nor_data <- ms_LFQ_nor_data[rownames(gp)]

ms_LFQ_nor_data <- na.omit(ms_LFQ_nor_data)


#################  30.	PCA and screeplot analysis
pca_scale <- PCAtools::pca(ms_LFQ_nor_data, scale = TRUE, metadata = gp)

pdf("PCA.pdf")
PCAtools::screeplot(pca_scale, components = getComponents(pca_scale)[1:10], axisLabSize = 10)
dev.off()


#################  31.	Generate a PCA plot
pdf(file ="PCA.pdf")
PCAtools::biplot(pca_scale, x="PC1", y="PC2",
                  colby = "condition", 
                  shape="condition",
                  shapekey = c(Pre_R=19, Pre_NR=20, Post_R = 12, Post_NR = 15, Normal = 18), 
                  pointSize = 2,
                  max.overlaps = 20,
                  #ellipse = TRUE, 
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  legendPosition = "right") 
dev.off()


#################  32.	Spearman correlation coefficient of plasma samples 
ms_LFQ_nor_data[is.na(ms_LFQ_nor_data)] <- 0
corvalue<-cor(ms_LFQ_nor_data, method = "spearman")

pdf(file="corrplot.pdf")
corrplot(corvalue, addgrid.col = "NA", addrect = 5, tl.cex = 0.5, type = "full", method = "color")
dev.off()


#################  33.	Pearson correlation coefficient of plasma proteins 
ms_LFQ_nor_data_pro <- t(ms_LFQ_nor_data)
ms_LFQ_nor_data_pro <- 2^ms_LFQ_nor_data_pro

corvalue_pro<-cor(ms_LFQ_nor_data_pro, method = "pearson")

pdf(file="corrplot_pro_heatmap.pdf")
pheatmap(corvalue_pro, fontsize = 6, fontsize_col = 1, fontsize_row = 1)
dev.off()


#################  34.	Consensus cluster analysis of all samples
ms_LFQ_nor_data <- as.matrix(ms_LFQ_nor_data)
concensus_results <- ConsensusClusterPlus(
                        ms_LFQ_nor_data,
                        maxK = 6,
                        reps = 50,
                        pItem = 0.8,
                        pFeature = 1,
                        clusterAlg = "hc",
                        distance = "spearman",
                        seed = 123,
                        title = "Consensus_All",
                        plot = "pdf",
                        writeTable = TRUE)
icl=calcICL(concensus_results, title = "Consensus_All", writeTable = TRUE, plot = "pdf")


#################  35.	Consensus cluster analysis of cancer samples
gp_cancer <- gp[gp$condition != "Normal", ]
ms_LFQ_nor_data_cancer <- ms_LFQ_nor_data[, rownames(gp_cancer)]
concensus_results <- ConsensusClusterPlus(
                        ms_LFQ_nor_data_cancer,
                        maxK = 6,
                        reps = 50,
                        pItem = 0.8,
                        pFeature = 1,
                        clusterAlg = "hc",
                        distance = "spearman",
                        seed = 123,
                        title = "Consensus_Cancer",
                        plot = "pdf",
                        writeTable = TRUE)
icl=calcICL(concensus_results, title = "Consensus_Cancer", writeTable = TRUE, plot = "pdf")


#################  36.	Read the survival information of the patients
survival_table <- xlsx::read.xlsx("survival.xlsx", header = TRUE, sheetIndex = 1)  
rownames(survival_table) <- survival_table[,1]
survival_table[,1] <- NULL
survival_table <- survival_table[order(rownames(survival_table)),]

head(rownames(survival_table))

head(survival_table)


#################  37.	Combine the protein expression matrix with the survival table
cancer51 <- ms_LFQ_nor_data[, rownames(survival_table)]  
head(rownames(cancer51))

colnames(cancer51)


## reverse the protein expression matrix
cancer51 <- t(cancer51)  

## combine the protein expression matrix with the survival table
cancer51 <- cbind(cancer51, survival_table)


#################  38.	Read the survival information of the patients
pro <- "FCN3_HUMAN"   ## The protein that will be analysed
pro_tab <- cbind(cancer51[,c("OS","OS_m","PFS","PFS_m")],
                  data.frame(Class = ifelse(cancer51[, pro] == "NA", "",
                     ifelse(cancer51[, pro] > median(cancer51[, pro], na.rm = T), 'High', "Low"))))
pro_sur_OS <- survival::survfit(Surv(OS_m, OS)~Class, data = pro_tab)
pro_sur_PFS <- survival::survfit(Surv(PFS_m, PFS)~Class, data = pro_tab)
sur_fun <- function(x){
  survminer::ggsurvplot(x, 
                        data = pro_tab,
                        legend.labs = levels(pro_tab[,'Class']),
                        legend.title = pro,
                        legend = 'top',
                        font.legend = c(12, "bold","black"),
                        fond.x = c(14, "bold", "black"),
                        fond.y = c(14, "bold", "black"),
                        surv.median.line = "hv",
                        conf.int = TRUE,
                        ggtheme = theme_bw(base_rect_size = 1, base_line_size = 0, 
                                           base_size = 15, base_family = 0),
                        pval = TRUE,
                        tables.y.text = FALSE, 
                        risk.table = TRUE,
                        risk.table.height = 0.27,
                        risk.table.col = "strata"
  )
}
pdf("survival.pdf")
sur_fun(pro_sur_OS)
sur_fun(pro_sur_PFS)
dev.off()


#################  END


