## Clear the environmental space here to start affresh 
rm(list = ls())

## Set the directory with the expression datasets here 
# Re-analysed samples with the new reference genome
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ11_star_output/star')

## Analysis (mapped on the pervious reference genome version)
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ10_star_output')

## This is the output folder for the final analysis
output_folder = '/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/finalAnalysis'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/Scripts/helper_function_and_env_setup.R")
source("/Users/gwk/Desktop/PhD/Scripts/star_output_to_count_matrix.R")

## Set up variables here 
## Read the mapped and quantified files into r for further analysis here 
file1 <- "1_S1_ReadsPerGene.out.tab"
dr <- getwd()

## load the data into R for further analysis here 
# This is output from star - these are files with ReadsPerGene.out.tab
countMatrix <- staroutput_preprocessing(dr, file1)

## Read in the data from the sample informatio here 
# If no sample information exists, you can create to match the countMatrix
samples <- read.csv("sample_information.csv", row.names = 1)

## Explore the loaded data here
head(countMatrix)
head(samples, 15)

## Fix the column names in the count matrix to match the sample information rownames
row_col_names <- as.vector(paste0("sample",rownames(samples)))
rownames(samples) <- row_col_names
colnames(countMatrix) <- row_col_names

## Makesure that sample rownames match the countmatrix colnames here
# The lines below must return true for all
# if false, please investigate the rownames match the countMatrix colnames
all(row.names(samples) %in% colnames(countMatrix))
all(row.names(samples) == colnames(countMatrix))

# Split the data so that mutant and wt 
# Mutant sample information
mut <- samples %>%
  filter(Genotype == "mutant") %>%
  dplyr::select(Group)

mutt <- c(3,4,8,9,10,11)
# Wild type sample information
wt <- samples %>%
  filter(Genotype == "wt") %>%
  dplyr::select(Group) 

## Split the matrix here to analyse each genotype seperately here
# Mutant count matrix
countMut <- countMatrix[,mutt]
names(countMut)

# WT countMatrix
countWT <- countMatrix[,-mutt]

# Check that all the rows and columns still match (sample rowname/count colnames)
all(row.names(mut) %in% colnames(countMut))
all(row.names(wt) %in% colnames(countWT))

## Change conditions to factor. This is the column or columns to be used by the
## model when buidling the object. Not Important, but if not done will return a
## warning message in the next step
mut$Group <- as.factor(mut$Group)
wt$Group <- as.factor(wt$Group)
samples$Group <- as.factor(samples$Group)

# Build a DESeq2 object here. The count data is given as a metrix, column names
# are the samples and the rownames are the Ensemble gene ids. This is important
# The design contains what makes up the model (negatve bionimial model in this case)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                              colData = samples,
                              design = ~ Group)

dds_mut <- DESeqDataSetFromMatrix(countData = countMut,
                                  colData = mut,
                                  design = ~ Group)

dds_wt <- DESeqDataSetFromMatrix(countData = as.matrix(countWT),
                                  colData = wt,
                                  design = ~ Group)
## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
pca_plot(dds_mut, mut)
pca_plot(dds_wt, wt)

## calculate the pca values here
vsd <- vst(dds)
plotPCA(vsd, intgroup="Group")

## Estimate the size factor here 
dds <- estimateSizeFactors(dds)
dds_mut <- estimateSizeFactors(dds_mut)
dds_wt <- estimateSizeFactors(dds_wt)
sizeFactors(dds_mut)
sizeFactors
normalised_all <- counts(dds, normalized = T)
normalised_counts_mut <- counts(dds_mut, normalized = T)
normalised_counts_wt <- counts(dds_wt, normalized = T)

## Export the data 
write.csv(normalised_all,
          file = paste0(output_folder,'/all_normalised.csv'))
write.csv(normalised_counts_mut,
          file = paste0(output_folder,'/mutant_normalised.csv'))
write.csv(normalised_counts_mut,
          file = paste0(output_folder,'/wt_normalised.csv'))

## Transform normalised data for visualisation here 
rld_mut <- rlog(dds_mut, blind = T)
rld_mut_mat <- assay(rld_mut)
pca_mut <- prcomp(t(rld_mut_mat))

## Wild type 
rld_wt <- rlog(dds_wt, blind = T)
rld_wt_mat <- assay(rld_wt)
pca_wt <- prcomp(t(rld_wt_mat))

## Mutants 
rld_mut <- rlog(dds_mut, blind = T)
rld_mut_mat <- assay(rld_mut)
pca_mut <- prcomp(t(rld_mut_mat))

# Create a data frame that can used going forward from here on 
count_mut_df <- cbind(mut,pca_mut$x)
count_wt_df <- cbind(wt, pca_wt$x)

## Visualise the data here 
ggplot(data = count_mut_df) +
  geom_point(aes(x=PC3, y=PC6, colour = Group), size=5) +
  theme_minimal() +
  labs(x = 'PC1: 68% variance',
       y = 'PC2: 22% varience') +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 15, vjust = 0.5),
    axis.text.y = element_text(angle = 90, vjust = 0.5),
    axis.title.y = element_text(size = 15, vjust = 0.5)
  )
  
## Perform hierarchical clustering 
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
rld_wt_corr <- cor(rld_wt_mat)

## Generat a heatmap plot here 
pheatmap(rld_cor)
pheatmap(rld_wt_corr)

## Perform DGE analysis here
dds_mut$Group <- relevel(dds_mut$Group, ref = "mut_unexposed")
dds_wt$Group <- relevel(dds_wt$Group, ref = "wt_unexposed")
dds$Group <- relevel(dds$Group, ref = "wt_unexposed")

dds_mut <- DESeq(dds_mut)
dds_wt <- DESeq(dds_wt)
dds <- DESeq(dds)

## Perform defferential gene expression analysis here
res_mut_treatment <- results(dds_mut, 
                         name =  "Group_mut_exposed_vs_mut_unexposed",
                         alpha = 0.05)

## Wild type 
res_wt_treatment <- results(dds_wt, 
                            name = "Group_wt_exposed_vs_wt_unexposed",
                            alpha = 0.05)

## All samples combined
res_all_mut_exposed <- results(dds,
                   name = "Group_mut_exposed_vs_wt_unexposed",
                   alpha = 0.05)

res_all_mut_unexposed <- results(dds,
                                  name = 'Group_mut_unexposed_vs_wt_unexposed',
                                  alpha = 0.05)
res_all_wt_exposed <- results(dds, 
                              name = 'Group_wt_exposed_vs_wt_unexposed',
                              alpha = 0.05)
## Check the summary of the DGE analysis output is here
summary(res_mut_treatment)
summary(res_wt_treatment)
summary(res_all_mut_exposed)
summary(res_all_mut_unexposed)
summary(res_all_wt_exposed)

## Perform shrunken analysis here 
res_mut_shrink <- lfcShrink(dds = dds_mut, 
                            coef = 2, 
                            res = res_mut_treatment)
summary(res_mut_shrink)

# Wild type animals 
res_wt_shrink <- lfcShrink(dds = dds_wt, 
                           coef = 2, 
                           res = res_wt_treatment)
summary(res_wt_shrink)

## All data combined here
res_all_mut_exposed_shrink <- lfcShrink(dds = dds,
                                        coef = 2,
                                        res = res_all_mut_exposed)

res_all_mut_unexposed_shrink <- lfcShrink(dds = dds,
                                          coef = 3,
                                          res = res_all_mut_unexposed)

res_all_wt_exposed_shrink <- lfcShrink(dds = dds,
                                       coef = 4,
                                       res = res_all_wt_exposed)

## Annotat the data before exporting the data 
res_all_wt_exposed_shrink$ENTREZID <- mapIds(org.Dr.eg.db,
                              keys = row.names(res_all_wt_exposed_shrink),
                              column = c("ENTREZID"),
                              keytype = "SYMBOL",
                              multiVals = "first")
  

## Wild type 
res_wt_shrink$ENTREZID <- mapIds(org.Dr.eg.db,
                                        keys = row.names(res_wt_shrink),
                                        column = c("ENTREZID"),
                                        keytype = "SYMBOL",
                                        multiVals = "first")

## Export data to a csv file for further analysis and inspection later
write.csv(res_all_wt_exposed_shrink,
          file = "/Users/gwk/Desktop/PhD/Data/PhD_data/Brain_Data_reanalysed/wt_exposed_lfc.csv")

write.csv(res_wt_shrink,
          file = "/Users/gwk/Desktop/PhD/Data/PhD_data/Brain_reanalysed/wildtype_reanalysed.csv")

################################################################################
############. Heatmap Plot of the DEG list here ################################
################################################################################
sig_gene <- function(datafframe, basemean = 100, lfc = 1.0) {
  sigs <- datafframe %>%
    as.data.frame() %>%
    filter(padj < 0.05)
  
  ## sig
  sig_df <- sigs[(sigs$baseMean > basemean) & abs((sigs$log2FoldChange >= lfc)),]
  sig_df <- na.omit(sig_df)
  
  ## Return value from the data frame here
  return(sig_df)
}

sig_df <- sig_gene(res_wt_shrink)

data_heatmap <- function(normalised_data, sample){
  sig_df <- sig_gene(res_wt_shrink)
  ## Return the z-score of the matrix here
  normalised_data <- normalised_data[rownames(sig_df),]
  mut_mat.z <- t(apply(normalised_data, 1, scale))
  colnames(mut_mat.z) <- sample$Group
  
  return(mut_mat.z)
  
}
data_heatmap(normalised_counts_wt, wt)
