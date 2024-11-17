## Clear the environment space here
rm(list = ls())

## Set the working directory 
setwd("/Users/gwk/Desktop/PhD/RNA_Sequencing/Bulk_RNA_Seq/PhD_RNA_Seq_Data/Re_run_counts/Final_Analysis")

## load data here
## Load the data to be further analysed

brain_mutant <- read.csv("Brain_Mutant_Exposed_vs_Mutant_Unexposed.csv", row.names = 1)
brain_wt <- read.csv("Brain_WT_Exposed_vs_WT_Unexposed.csv", row.names = 1)

## Set the working directory here 
setwd("/Users/gwk/Desktop/PhD/Data/PhD_data/")

# Check the working directory
dir()

## Set up variables 
file1 <- "1_S1_ReadsPerGene.out.tab"
dr <- getwd()

## load the data into R for further analysis here 
countMatrix <- staroutput_preprocessing(dr, file1)

# Edit the names of the col in the data
# colnames(countMatrix) <- sub("*_S", "",colnames(countMatrix))

## Read in the data from the sample informatio here 
samples <- read.csv("sample_information.csv", row.names = 1)

## fix the column and row names in the sample and count matrix 
row_col_names <- as.vector(paste0("sample",rownames(samples)))
rownames(samples) <- row_col_names
colnames(countMatrix) <- row_col_names

## Confirm that rownames in the samples match the colnames in the count matrix
all(rownames(samples)) %in% all(colnames(countMatrix))

# Split the data so that mutant and wt 
mut <- samples %>%
  filter(Genotype == "mutant")
wt <- samples %>%
  filter(Genotype == "wt")

mut <- mut[c(1,3,5,2,4,6),]
wt <- wt[c(1,3,6,5,2,4),]

## Final count matrix for Differential gene expression analysis
mut_norm <- countMatrix %>%
  dplyr::select(rownames(mut))

wt_norm <- countMatrix %>%
  dplyr::select(rownames(wt))

## Build the object here
# dds for the mutant groups
dds_mut <-  DESeqDataSetFromMatrix(countData = as.matrix(mut_norm),
                                   colData = mut,
                                   design = ~ Group)

# dds fores# dds for the wild type group
dds_wt <- DESeqDataSetFromMatrix(countData = as.matrix(wt_norm),
                                 colData = wt,
                                 design = ~ Group)

dds_all <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                                  colData = sample_information,
                                  design = ~ Group)

## Perform differential gene analysis here
dds <- DESeq(dds_mut)
res <- results(dds, alpha = 0.05)
summary(res)
res[res$pvalue < 0.05,]
plot(sizeFactors(dds_mut), colSums(counts(dds_mut)))

## Only mutant heatmap
sig_mut_df <- df[(df$baseMean > 100) & abs(df$log2FoldChange > 1.2),]
sig_mut_df <- na.omit(sig_mut_df)

dds_mut <- estimateSizeFactors(dds_mut)
mat_mut <- counts(dds_mut, normalized =T)[sig_mut_df,]
mat_mut <- t(apply(mat_mut, 1, scale))
colnames(mat_mutt) <- mut$Group

## Revove any NaN from the matrix here
mat_mut <- mat_wt[-c(33),]
sig_wt_df <- sig_wt_df[-c(33),]

## Plot the heat map here 

## Only wt heatmap here
df <- res_mut_treatment
sig_wt_df <- df[(df$baseMean > 100) & abs(df$log2FoldChange > 1.2),]
sig_wt_df <- na.omit(sig_wt_df)
wt_norm <- estimateSizeFactors(dds_wt)
mat_wt <- counts(wt_norm, normalized = T)[rownames(sig_wt_df),]
mat_wt <- t(apply(mat_wt, 1, scale))
colnames(mat_wt) <- wt$Group
mat_wt <- mat_wt[-c(33),]
sig_wt_df <- sig_wt_df[-c(33),]

## Heatmap 
figure_heatmap <- function(data_frame, sigs){
  plot1 <- Heatmap(data_frame, name = "Z - Score", row_km = 2, column_km = 2,row_labels = rownames(sigs),
                   column_labels = colnames(data_frame),
                   top_annotation = HeatmapAnnotation(data_frame = 1:dim(data_frame)[2], bar1 = anno_points(runif(dim(data_frame)[2]))),
                   right_annotation = rowAnnotation(data_frame = dim(data_frame)[1]:1, 
                                                    bar2 = anno_barplot(runif(dim(data_frame)[1])))
  )
  
  return(plot1)
}

figure_heatmap(normalised_counts_wt, sig_et_df)

### Gene ontology analysis output and visualisation 
## Read the data into R here 
setwd("/Users/gwk/Desktop/PhD/Data/PhD_data/analysis/")
mut_go <- read.csv("brain_mutant_gse_GO_simplified.csv")
wt_go <- read.csv("brain_wt_gse_simplified.csv")

## Explore the data here 
head(mut_go)
names(mut_go)

## Wild type 
wt_go$Description[grep("synap", wt_go$Description)]
wt_go$Description[grep("calcium", wt_go$Description)]
wt_go$Description[grep("oxi", wt_go$Description)]

## Mutant Terms 
gse_simplified$Description[grep("synap", gse_simplified$Description)]
gse_simplified$Description[grep("calcium", gse_simplified$Description)]
gse_simplified$Description[grep("oxi", gse_simplified$Description)]

## 
s_transmission <- mut_go[grep("calcium-ion regulated exocytosis", mut_go$Description),]
s_transmission |>
  dplyr::select(ONTOLOGY, ID, Description, NES, qvalue) |>
  tibble::tibble() |>
  print(n = 25)


wt_go |>
  filter(Description == wt_go$Description[grep("synap", wt_go$Description)])


syn <- as.data.frame(g_name)

syn %>%
  arrange(desc(g_name))

## Export the lists here 
write.csv(syn, 
          file = "calmodulin binding.csv")

## Visualise GO Terms
dff2 <- mut_go |>
  as.data.frame() |>
  filter(ONTOLOGY == 'BP') |>
  filter(NES < -1.9) |>
  # Split the leading edge into its components for visualisation here
  separate( 
           leading_edge, 
           into = c('Tags', "lists","signal"),
           sep = ",") |>
  mutate(Tags = gsub("tags=", "", Tags),
         GeneRatio = as.numeric(gsub("%","", Tags))/100) |>
  arrange(desc(NES)) |>
  dplyr::select(ID, Description, enrichmentScore, NES, qvalue, setSize, GeneRatio) 
  attach(ndff)

  ggplot(data = ndff, aes(x = NES, y = Description)) +
  geom_segment(x = NES, y = Description, xend = NES, y = 0, yend = Description) +
  geom_point(aes(size = GeneRatio, colour = )) +
  labs(
    x = "Normalised Enrichment Score",
    y = "Gene ontology terms"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(size = 12.5, face = 'bold'),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    plot.caption = element_text(hjust=0.02, size=rel(2.2))
  ) 
dff2
ndff <- rbind(dff2, dff1)

## Tables of data
# Top 10 upregulated terms
top10_up <- mut_go |>
  filter(NES > 1.5 | NES < -1.5) |>
  # Split the leading edge into its components for visualisation here
  separate( 
    leading_edge, 
    into = c('Tags', "lists","signal"),
    sep = ",") |>
  mutate(Tags = gsub("tags=", "", Tags),
         `Gene Ratio` = as.numeric(gsub("%","", Tags))/100) |>
  arrange(desc(NES)) |>
  dplyr::select(Description, enrichmentScore, NES, p.adjust, setSize, `Gene Ratio`) |>
  head(10) |> tibble()

# Top 10 downregulated terms
top10_down <- mut_go |>
  filter(NES > 1.5 | NES < -1.5) |>
  # Split the leading edge into its components for visualisation here
  separate( 
    leading_edge, 
    into = c('Tags', "lists","signal"),
    sep = ",") |>
  mutate(Tags = gsub("tags=", "", Tags),
         `Gene Ratio` = as.numeric(gsub("%","", Tags))/100) |>
  arrange(desc(NES)) |>
  dplyr::select(Description, enrichmentScore, NES, p.adjust, setSize,`Gene Ratio`) |>
  tail(10) |>
  tibble()

## Export the table for insertion into the presentation
write.csv(top10_up,
          file = "mut_upregulated_terms.csv")

write.csv(top10_down,
          file = "mut_downregulated_terms.csv")

## Check the score
mutate(mut_go, qscore = -log(p.adjust, base = 10)) |>
  barplot(x = "qscore")

## KEGG visualisation and analysis here 
# Read the data into here 
kegg_pathways <- read.csv()

















