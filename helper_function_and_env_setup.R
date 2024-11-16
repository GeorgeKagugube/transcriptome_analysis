## Function for preparing gene lists for gsea or over representation analysis
## Developed by George William Kagugube 
## Date: 29 Jan 2024 
## Usage, read teh file into R preferably using source (file_name)
## prepare an object from DESeq of the differentially expressed genes
## This is given as an  data frame to the function here

# Load the libraries that will be needed for the analysis herein 
requiredPackages <- c('ggvenn', 'dplyr', 'tidyr', 'data.table', 'ggplot2',
                      'pheatmap', 'RColorBrewer', 'ComplexHeatmap', 'circlize',
                      'EnhancedVolcano', 'gridExtra', 'grid', 'enrichplot',
                      'clusterProfiler', 'org.Dr.eg.db', 'GOSemSim','DOSE',
                      'ggupset','AnnotationDbi', 'AnnotationHub', 'DESeq2', 
                      "pathview","ReactomePA", "gtExtras")

##################### Loading the required libraries here ######################
lapply(requiredPackages, library, character.only=TRUE)

################################################################################
################### Functions start from here. #################################
################################################################################
transcriptome_annotation <- function(df){
  # Inputs: a dataframe from DESeqs (containing ENSEMBL IDs and log fold changes)
  # Inputs: choose an organism to work with, default is danio - Zebrafish
  # Outputs: annotated dataframe with gene symbols and entrezids
  # The Ensembl IDs need to be the row names for this function to work
    print("Annotating data with Zebrafish gene symbols and EntrezIDs......")
    df$Symbols <- mapIds(org.Dr.eg.db,
                         keys = rownames(df),
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
    
    df$Entrezid <- mapIds(org.Dr.eg.db,
                          keys = rownames(df),
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")
    
    ## Return the annotated data object here 
    return(df)
}

ora_gene_list <- function(df, cutoff=1.0){
  # Function that creates a list of genes to use in functional annotation
  # This comuptes over-representation of genes from our list into 
  # Inputs: A dataframe from DESeq, with the rows names == Ensembl gene IDs
  # A cuttoff value for logfoldchange, default is 1.0
  # Output1: Gene Universe (aka reference/Background) - all the genes identified 
  # in the analysis
  # Output2: A vector of gene list with significant expression and lfc > cuttoff
  
  # Extract the log2foldchange from the lfcshrunk list
  gene_universe <- df$log2FoldChange
  
  ## Name the vector
  names(gene_universe) <- rownames(df)
  
  # Remove any na values from the list
  gene_list <- na.omit(gene_universe)
  
  # Sort the list here 
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Extract significant genes here
  sig_gene_df = subset(df, padj < 0.05)
  
  # Extract the log2fold change here
  genes <- sig_gene_df$log2FoldChange
  
  # Name the vector here
  names(genes) <- rownames(sig_gene_df)
  
  ## Remove any missing log2fold values here
  genes <- na.omit(genes)
  
  # filter by log2fold value here
  genes <- names(genes)[abs(genes) > cutoff]
  
  genes_ids <- list(gene_list, genes)
  
  return (genes_ids)
}

ora_gene_list_entrezid <- function(df, cutoff=1.0){
  # Function that creates a list of genes to use in functional annotation
  # This comuptes over-representation of genes from our list into 
  # Inputs: A dataframe from DESeq, with the rows names == Ensembl gene IDs
  # A cuttoff value for logfoldchange, default is 1.0
  # Output1: Gene Universe (aka reference/Background) - all the genes identified 
  # in the analysis
  # Output2: A vector of gene list with significant expression and lfc > cuttoff
  
  # Extract the log2foldchange from the lfcshrunk list
  gene_universe <- df$log2FoldChange
  
  ## Name the vector
  names(gene_universe) <- df$ENTREZID
  
  # Remove any na values from the list
  gene_list <- na.omit(gene_universe)
  
  # Sort the list here 
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Extract significant genes here
  sig_gene_df = subset(df, padj < 0.05)
  
  # Extract the log2fold change here
  genes <- sig_gene_df$log2FoldChange
  
  # Name the vector here
  names(genes) <- sig_gene_df$ENTREZID
  
  # filter by log2fold value here
  genes <- names(genes)[abs(genes) > cutoff]
  
  ## Remove any missing log2fold values here
  genes <- na.omit(genes)
  
  genes_ids <- list(gene_list, genes)
  
  return (genes_ids)
}

gsea_gene_generater_func <- function(df){
  # Function that computes a list of genes to use in functional annotation
  # Inputs: A dataframe from DESeq, with the rows names == Ensembl gene IDs
  # Outputs: A ranked vector of gene list from high to low 
  gsea_gene_list <- df$log2FoldChange
  names(gsea_gene_list) <- df$Entrezid
  gsea_list <- na.omit(gsea_gene_list)
  gsea_list <- sort(gsea_list, decreasing = TRUE)
  
  return(gsea_list)
}

reactome_gsea_gene_generater_func <- function(df){
  # Function that computes a list of genes to use in functional annotation
  # Inputs: A dataframe from DESeq, with the rows names == Ensembl gene IDs
  # Outputs: A ranked vector of gene list from high to low 
  gsea_gene_list <- df$log2FoldChange
  names(gsea_gene_list) <- df$Entrezid
  gsea_list <- na.omit(gsea_gene_list)
  gsea_list <- sort(gsea_list, decreasing = TRUE)
  
  return(gsea_list)
}

gene_names <- function(df, term){
  # The function takes in a GO term dataframe from clusterprofiler
  # Inputs: 
          # 1. Dataframe containing gseGO output
          # 2. term is any term that one needs to 
  # 
  ensembl_id_list <- unlist(strsplit(gsub("/",",",df[df$Description == term,]$core_enrichment),","))
  gene_symbols <- mapIds(org.Dr.eg.db,
                   keys = ensembl_id_list,
                   column = c("SYMBOL"),
                   keytype = c("ENSEMBL"),
                   multiVals = "first")
  
  return(gene_symbols)
}

ora_kegg_gene_list <- function(df, cuttoff = 1.0){
  ## ORA KEGG Pathway
  ## Create a vector of the gene inverse
  kegg_gene_list <- df$log2FoldChange
  
  ## Name the vector with ENTREZIDS
  names(kegg_gene_list) <- df$Entrezid
  
  # Remove any NA values
  kegg_gene_list <- na.omit(kegg_gene_list)
  
  # Sort the lost here 
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
  
  # Extract significant results from df (create the marginal gene list)
  kegg_sig_genes_df <- subset(df, padj < 0.05)
  
  ## From the df above, extract the gene list
  kegg_genes <- kegg_sig_genes_df$log2FoldChange
  
  # name the vector above here
  names(kegg_genes) <- kegg_sig_genes_df$Entrezid
  
  ## filter by logfoldchange greater than 1 here
  kegg_genes <- names(kegg_genes)[abs(kegg_genes) > cuttoff]
  
  # Remove all the NA values from the vector
  kegg_genes <- na.omit(kegg_genes)
  
  return(kegg_genes)
}

gseGO_gene_name_decoding <- function(df, word = "synap"){
  
  # Return the gene names
  gene_symbols <- df$Description[grep(word, df$Description)]
  
  # Return the genes
  return(gene_symbols)
}

select_genes_of_go_term <- function(dataframe,esembl_ids){
  # This function allows the user to exctract gene names and logFoldChange
  # Inputs:
  # 1. a dataframe, whose rownames are esemble gene ids 
  # 2. a list of ensemble gene ids
  
  # Outputs:
  # a table of genes names and their corresponidng logFoldChange values
  
  data_out <- dataframe[esembl_ids,] %>%
    dplyr::select(Symbols,log2FoldChange)
  
  return(data_out)
}

sample_dist <- function(df, sample_id){
  
  # To plot the distribution of the data of any chosen sample
  # Inputs: 
          # 1. df, a normalised data frame containing the read count data
          # 2. sample_id, any chosen column of the count data frame here
  # Outputs: Visualisation of the sample distribution 
  plot1 <- ggplot(data = df) +
    geom_histogram(aes(x=df$sample_id), stat = "bin", bins = 200) +
    # zoom in to see the distribution more clearly
    xlim(-2, 50) +
    xlab("Raw expression counts") +
    ylab("Number of genes") +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 30),
      axis.text = element_text(family = "Arial" ,face = "bold",size = 30),
      axis.text.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.x = element_text(family = "Arial", face = "bold", size = 30)
    )
  
  return(plot1)
}

mean_var_plot <- function(df){
  mean_counts <- apply(countData[, 1:length(colnames(df))], 1, mean)
  variance_counts <- apply(df[, 1:length(colnames(df))], 1, var)
  df_count <- data.frame(mean_counts, variance_counts)
  
  ## Plot here
  plot2 <- ggplot(df_count) +
    geom_point(aes(x=mean_counts, y=variance_counts)) + 
    geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
    scale_y_log10() +
    scale_x_log10() +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 30),
      axis.text = element_text(family = "Arial" ,face = "bold",size = 30),
      axis.text.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.x = element_text(family = "Arial", face = "bold", size = 30)
    )
  return(plot2)
}

pca_plot <- function(dds, samples){
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  # Extract the rlog matrix from the object
  rld_mat <- assay(rld)
  colnames(rld_mat) <- samples$Group
  pca <- prcomp(t(rld_mat))
  
  z <- plotPCA(rld, intgroup="Group")
  z + geom_label(aes(label = samples)) +
    theme_classic() +
    theme(
      #text = element_text(family = "Arial", face = "bold", size = 30),
      #axis.text = element_text(family = "Arial" ,face = "bold",size = 30),
      axis.text.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.text.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.x = element_text(family = "Arial", face = "bold", size = 30)
    )
  
  return(z)
}

################################################
############# Gene Counts plots   #############
################################################
gene_count <- function(df, gene) {
  plot_gene <- df %>%
    filter(Name == gene) %>%
    dplyr::select(starts_with("Sa")) %>%
    gather(key = "Samples", value = "Counts")  %>%
    mutate(Treatment = c("wt_unexposed","wt_unexposed","wt_unexposed",
                         "mut_unexposed","mut_unexposed","mut_unexposed",
                         "wt_exposed","wt_exposed","wt_exposed",
                         "mut_exposed","mut_exposed","mut_exposed")) %>%
    ggplot(aes(x = Treatment, y = Counts)) +
    geom_jitter(position = position_dodge(0.2)) +
    geom_boxplot(alpha=0.2) +
    scale_x_discrete(limits = c("wt_unexposed","wt_unexposed","wt_unexposed",
                                "wt_exposed","wt_exposed","wt_exposed",
                                "mut_unexposed","mut_unexposed","mut_unexposed",
                                "mut_exposed","mut_exposed","mut_exposed")) +
    labs(
      title = gene,
      y = "Normalised Counts",
      x = "Treatment (µM)"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(family = "Times-Roman", size = 15),
      axis.text.x = element_text(family = "Times-Roman",size = 15),
      axis.title.y = element_text(family = "Times-Roman",size = 20),
      axis.title.x = element_text(family = "Times-Roman",size = 20),
      legend.text = element_text(family = "Times-Roman",size = 15),
      legend.title = element_text(family = "Times-Roman",size = 20)
    )
  
  return(plot_gene)
}

####### valcano plots
volcanoPlot <- function(df){
  plot_1 <- EnhancedVolcano(df,
                            lab = rownames(df),
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            axisLabSize = 30,
                            xlim = c(-3, 3),
                            ylim = c(0, 30),
                            pointSize = 2.0,
                            labSize = 8.5,
                            labFace = "plain",
                            title = 'DESeq2 results',
                            subtitle = 'Differential expression',
                            caption = 'FC cutoff, 1.333; p-value cutoff, 0.05',
                            legendPosition = "right",
                            legendLabSize = 30,
                            pCutoff = 0.05,
                            FCcutoff = 1.0)
  
  return(plot_1)
}

venPlot <- function(dgeset){
  p3 <- ggvenn(dgeset, 
               set_name_size = 14, 
               text_size = 7,
               fill_alpha = 0.25,
               stroke_size = 1.5,
               show_outside = "auto")
  
  return(p3)
}

# Function for generating the 
enrich_ont <- function(dge_file, method="gsea"){
  
  if (method == "gsea"){
    # Display a message of what analysis is being conducted here
    print ("You're analysing your gene list using gsea")
    # Select the fold change from the 
    original_gene_list <- dge_file$log2FoldChange
    names(original_gene_list) <- rownames(dge_file)
    original_gene_list <- na.omit(original_gene_list)
    gene_list <- sort(original_gene_list, decreasing = TRUE)
    
    return(gene_list)
  }else{
    print ("You're analysing your gene list using over representation")
    # we want the log2 fold change 
    
    gene_list <- dge_file[(abs(dge_file$log2FoldChange > 0.5) & dge_file$padj < 0.05),]
    genes <- rownames(gene_list)
    
    return(genes)
  }
}

### Pathway analysis gene list preperation 
kegg_pathway_analysis <- function(dge_file, method = "gsea"){
  
  # Select the fold change from the 
  original_gene_list <- dge_file$log2FoldChange
  names(original_gene_list) <- rownames(dge_file)
  
  if (method == "gsea"){
    print ("You're analysing your gene list using gsea")
    ## Convert the ensembl IDs to kegg API naming convention 
    ids <- bitr(geneID = names(original_gene_list),
                     fromType = "ENSEMBL",
                     toType = "ENTREZID",
                     OrgDb = "org.Dr.eg.db",
                     drop = TRUE)
    
    # Remove any duplicates from the list above 
    dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),]
    
    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    df2 = dge_file[rownames(dge_file) %in% dedup_ids$ENSEMBL,]
    
    # Create a new column in df2 with the corresponding ENTREZ IDs
    rownames(dedup_ids) <- dedup_ids$ENSEMBL 
    
    df2 <- merge(dedup_ids, dge_file, by="row.names")
    
    ## Select the desired genes here 
    kegg_gene_list <- df2$log2FoldChange
    names(kegg_gene_list) <- df2$ENTREZID
    kegg_gene_list <- na.omit(kegg_gene_list)
    kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
    
    kegg_gene_list <- kegg_gene_list[!duplicated(kegg_gene_list)]
    
    return(kegg_gene_list)
    
  }else{
    print ("You're analysing your gene list using over representation")
    ids<-bitr(names(original_gene_list), 
              fromType = "ENSEMBL", 
              toType = "ENTREZID", 
              OrgDb="org.Dr.eg.db")
    
    dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
    df2 = dge_file[rownames(dge_file) %in% dedup_ids$ENSEMBL,]
    # Create a new column in df2 with the corresponding ENTREZ IDs
    df2$ENTREZID = dedup_ids$ENTREZID
    
    # Create a vector of the gene unuiverse
    kegg_gene_list <- df2$log2FoldChange
    
    # Name vector with ENTREZ ids
    names(kegg_gene_list) <- df2$ENTREZID
    
    # omit any NA values 
    kegg_gene_list<-na.omit(kegg_gene_list)
    
    # sort the list in decreasing order (required for clusterProfiler)
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
    
    # Exctract significant results from df2
    kegg_sig_genes_df = subset(df2, padj < 0.05)
    
    # From significant results, we want to filter on log2fold change
    kegg_genes <- kegg_sig_genes_df$log2FoldChange
    
    # Name the vector with the CONVERTED ID!
    names(kegg_genes) <- kegg_sig_genes_df$ENTREZID
    
    # omit NA values
    kegg_genes <- na.omit(kegg_genes)
    
    # filter on log2fold change (PARAMETER)
    kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 0.5]
    
    return(kegg_genes)
  }
}

##### Functional analysis by calling gse
gseFunc <- function(gene_list, item="ALL"){
  gse <- gseGO(geneList = gene_list,
               ont = item,
               keyType = "ENSEMBL",
               minGSSize = 5,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               eps = 0,
               OrgDb = "org.Dr.eg.db",
               pAdjustMethod = "fdr")
  return(gse)
}

oraFunc <- function (gene_list, item = "ALL") {
  go_enrich <- enrichGO(gene = gene_list,
                        OrgDb = "org.Dr.eg.db", 
                        keyType = 'ENSEMBL',
                        ont = item,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = 0.10)
  return(go_enrich)
}

pathway_ont <- function(genelist, mode = "gse"){
  if (mode == "gse"){
    print("Using the entire gene list for analysis here")
    kk <- gseKEGG(geneList     = genelist,
                  organism     = "dre",
                  minGSSize    = 3,
                  maxGSSize    = 800,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr",
                  keyType       = "ncbi-geneid")
  }else{
    kk <- enrichKEGG(gene=genelist, 
                     organism="dre", 
                     pvalueCutoff = 0.05, 
                     keyType = "ncbi-geneid")
  }
  return(kk)
}

### Visualisation functions
plotting_GO <- function(obj, mode = "gsea") {
  
  if (mode == "gsea"){
    plot <- dotplot(obj, 
                    color = "p.adjust",
                    showCategory=10, 
                    font.size = 20,
                    split=".sign") + facet_grid(.~.sign)
    
    return(plot)
  }else{
    plot <- dotplot(obj, 
                    color = "p.adjust",
                    showCategory=10, 
                    font.size = 20)
    
    return(plot)
  }
  plot <- dotplot(obj, 
                  color = "p.adjust",
                  showCategory=10, 
                  font.size = 17,
                  split=".sign") + facet_grid(.~.sign)
  
  return(plot)
}

emap_plotting <- function(obj){
  ego2 <- pairwise_termsim(obj,
                           method = "JC")
  plot <- emapplot(ego2)
  
  return(plot)
}
