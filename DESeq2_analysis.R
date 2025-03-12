#Set working directory
setwd("//wsl.localhost/Ubuntu/home/jyothi/Github_projects/RA_Mar5")

#Install required packages
BiocManager::install("GEOquery", force = TRUE)
library(GEOquery)
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("limma", force = TRUE)
library(limma)
BiocManager::install("pheatmap")
library(pheatmap)
library(ggplot2)
#Load the GEO datasets
#GSE82107 (OA dataset)
gse82107 <- getGEO("GSE82107", GSEMatrix = TRUE)
#GSE77298 (RA dataset)
gse77298 <- getGEO("GSE77298", GSEMatrix = TRUE)
#GSE55235 (contains both OA and RA samples)
gse55235 <- getGEO("GSE55235", GSEMatrix = TRUE)
#GSE1919 (contains both OA and RA samples)
gse1919 <- getGEO("GSE1919", GSEMatrix = TRUE)

#Extract the expression data
exprs_gse82107 <- exprs(gse82107[[1]])
exprs_gse77298 <- exprs(gse77298[[1]])
exprs_gse55235 <- exprs(gse55235[[1]])
exprs_gse1919 <- exprs(gse1919[[1]])

#Extract the platform data
gpl_82107 <- Table(getGEO("GPL570"))
gpl_77298 <- Table(getGEO("GPL570"))
gpl_55235 <- Table(getGEO("GPL96"))
gpl_1919 <- Table(getGEO("GPL91"))

#Combine the expression & platform data into a list
expression_matrices <- list(exprs_gse82107, exprs_gse77298, exprs_gse55235, exprs_gse1919)
platform_files <- list(gpl_82107, gpl_77298, gpl_55235, gpl_1919)
processed_exp_list <- list()

#Loop through each expression file and its corresponding platform file
for (i in 1:length(expression_matrices)) 
{
  #Extract current expression matrix and platform data
  exp <- expression_matrices[[i]]
  plate <- platform_files[[i]]
  #Extract the ID and GENE SYMBOL columns from the platform file
  ID <- data.frame(ID_REF = plate$ID, Gene_Symbol = plate$`Gene Symbol`)
  #Remove the trailing spaces
  ID$Gene_Symbol <- data.frame(sapply(ID$Gene_Symbol, function(x) unlist(strsplit(x, "///"))[1]), stringsAsFactors=F)[,1]
  #Store it as a data frame
  exp <- as.data.frame(exp)
  #Make sure the rownames match the column names
  exp<-tibble::rownames_to_column(exp, "ID_REF")
  #Merge both the files
  exp <- merge(exp, ID, by='ID_REF')
  #Remove leading and trailing whitespace from data
  exp[, grep("Gene_Symbol", colnames(exp))] <-  trimws(exp[, grep("Gene_Symbol", colnames(exp))])
  #Replace empty values with NA
  exp[exp == ""] <- NA
  #Remove data with missing GENE_SYMBOL
  exp <- na.omit(exp)  
  #Save the processed data into a list
  processed_exp_list[[i]] <- exp
  #Save it as .csv file
  write.csv(exp, paste0("Processed_Expression_Dataset_", i, ".csv"), row.names = FALSE)
}

#Define file names for the count matrices and metadata
expression_files <- list("Processed_Expression_Dataset_1.csv", 
                         "Processed_Expression_Dataset_2.csv", 
                         "Processed_Expression_Dataset_3.csv", 
                         "Processed_Expression_Dataset_4.csv")
metadata_files <- list("Metadata_Dataset_1.csv", 
                       "Metadata_Dataset_2.csv", 
                       "Metadata_Dataset_3.csv", 
                       "Metadata_Dataset_4.csv")

#Open a PDF file to save all MA and volcano plots
pdf("All_MA_Volcano_Heatmap_plots.pdf")

#Loop through each dataset and metadata pair
for (i in 1:length(expression_files)) 
{
  #Read the expression dataset
  ped <- read.csv(expression_files[[i]])
  #Make duplicate rows unique by gene symbol
  rownames(ped) <- make.unique(as.character(ped$Gene_Symbol))
  #Remove the first and last column to get the count matrix
  count_matrix <- ped[, 2:(ncol(ped) - 1)]
  #Read the metadata file
  metadata <- read.csv(metadata_files[[i]])
  #Create the DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = metadata, design = ~ title)
  #Run DESeq2 analysis
  dds <- DESeq(dds)
  #Extract results
  res <- results(dds)
  #Convert the results to a data frame
  res_df <- as.data.frame(res)
  #Add a column to classify significant vs non-significant
  res_df$sig <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, "Significant", "Not Significant")
  #Filter for only significant genes
  res_sig_df <- res_df[res_df$sig == "Significant", ]
  #Save the significant DEGs to a CSV file
  write.csv(as.data.frame(res_sig_df), file = paste0("Significant_DEGs", i, ".csv"), row.names = TRUE)
  #Plot and add the MA plot to the PDF
  par(mfrow=c(1,1))  # Ensure plots are printed one per page
  DESeq2::plotMA(res, ylim = c(-5, 5), main = paste("MA Plot - Dataset", i))
  #Generate the volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(col = sig)) +
    scale_color_manual(values = c("red", "black")) +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = paste("Volcano Plot - Dataset", i))
  #Print the volcano plot to the PDF
  print(volcano_plot)
}
#Close the PDF device to save all plots
dev.off()


#Define file names for the processed expression data and updated DEGs
exp_files_with_symbols <- c("Gene_symbol_Processed_Expression_Dataset_1.csv", 
                            "Gene_symbol_Processed_Expression_Dataset_2.csv", 
                            "Gene_symbol_Processed_Expression_Dataset_3.csv", 
                            "Gene_symbol_Processed_Expression_Dataset_4.csv")

updated_DEGs <- c("Updated_Significant_DEGs1.csv",
                  "Updated_Significant_DEGs2.csv", 
                  "Updated_Significant_DEGs3.csv", 
                  "Updated_Significant_DEGs4.csv")

#Open a PDF file to save all heatmaps
pdf("All_Heatmaps.pdf", width = 10, height = 8)
#Loop through the files
for (i in 1:4) 
{
  #Read the expression files
  u1 <- read.csv(exp_files_with_symbols[i], row.names = 1)
  #Read the DEG file
  u2 <- read.csv(updated_DEGs[i], row.names = 1)
  #Convert the rownames as a separate column - in this case, gene symbols
  u2 <- cbind(Gene_Symbol = rownames(u2), u2)
  #Remove the rownames
  rownames(u2) <- NULL
  #Find common genes between both the files
  common_genes <- intersect(u1$Gene_Symbol, u2$Gene_Symbol)
  #Print the number of common genes
  print(paste("Number of matching genes:", length(common_genes)))
  #Subset the genes
  expr_deg <- u1[u1$Gene_Symbol %in% u2$Gene_Symbol, ]
  #Print the subset genes
  print(head(expr_deg))
  #Take average if duplicate gene entries are present by grouping all gene symbols and averaging the numeric columns
  expr_deg_avg <- expr_deg %>%
    group_by(Gene_Symbol) %>%                    
    summarise(across(where(is.numeric), mean))
  #Convert to a matrix
  expr_matrix <- as.matrix(expr_deg_avg[, -1])  # Remove Gene_Symbol column
  #Perform scaling
  expr_scaled <- t(scale(t(expr_matrix)))
  #Convert back to a data frame
  expr_scaled <- as.data.frame(expr_scaled)
  #Generate heatmap
  pheatmap(expr_scaled,
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = TRUE, show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           fontsize_row = 6,
           main = paste("Heatmap for ", i))
}
#Close the PDF file
dev.off()
