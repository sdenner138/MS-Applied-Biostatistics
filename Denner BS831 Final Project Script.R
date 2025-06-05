## BS831 Final Project
## Sara (Shiv) Denner
## Due 7 May 2025


library(Biobase)
library(dplyr)

setwd('C:/Users/sarad/OneDrive/Desktop/School/MSAB/Spring 2025/BS831/Final Project')

# read in data
fp <- readRDS("HNSC_htseq_raw_counts_AEvsG1vsG3.rds")
class(fp)
dim(fp)

#34422 genes, 120 samples


################## Pre-processing & QC

##### subset to keep only samples with grade g1 or g3


keep <- pData(fp)$grade == "g3" | pData(fp)$grade == "g1"

fp_filtered <- fp[,keep]

table(fp_filtered$grade)

###### remove genes with total counts equal to zero

n_zero <- rowSums(exprs(fp_filtered)) == 0
table(n_zero)


fp_filtered <- fp_filtered[!n_zero,]


# find distribution of expression values
hist(exprs(fp_filtered))
dim(fp_filtered)
# After filtering out zero counts and reducing to g3 and g1,
# 32208 genes and 80 samples.

############################ Filtering out Low Counts & Running DESeq2

library(DESeq2)

# format dataset for DESeq2

col_data_fp <- data.frame(
  condition = factor(as.character(pData(fp_filtered)[, "grade"]), levels = c("g1", "g3"))
)

dds_fp <- DESeqDataSetFromMatrix(
  countData = exprs(fp_filtered), 
  colData = col_data_fp, 
  design = ~ condition
)

## according to bioconductor vignette, pre-filter after creating dds object.
## Perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples
## since the smallest number of samples is 40 after filtering out the AE condition, we will
## ensure that rows have a count of at least 10 for 40 samples.

smallestGroupSize <- 40
keep_dds <- rowSums(counts(dds_fp) >= 10) >= smallestGroupSize
dds_fp <- dds_fp[keep_dds, ]

# check how many genes remain after pre-filtering
nrow(dds_fp)

# after filtering out low counts, 16,967 genes remain.


# run DESeq2
dds_res_fp <- DESeq(dds_fp)

# extract results
res_deseq_fp <- results(dds_res_fp)

head(res_deseq_fp)

# number of significant genes
table(res_deseq_fp$padj < 0.05)
# 4007 genes are significant

res_deseq_sort_fp <- res_deseq_fp[order(res_deseq_fp$padj), ]
head(res_deseq_sort_fp)

### Top 3 Most Highly Significantly Differentially Expressed genes
### gene IDs: ENSG00000197915, ENSG00000134765, ENSG00000176075
### raw p-values: 4.34947e-27, 1.17913e-25, 1.11990e-21
### FDR adjusted p values: 7.37975e-23, 1.00032e-21, 6.33380e-18
### log2FoldChange: -5.33101, -5.44052, -5.13362, means the top differentially expressed genes are more highly expressed in G1 than in G3
### large positive stats mean the gene is upregulated in G3, large negative stats mean gene is upregulated in G1


# gene symbols for top 3 genes above
# HRNR DSC1 LINC00302
fData(fp_filtered)[rownames(res_deseq_sort_fp)[1:3], ]

#MA Plot
DESeq2::plotMA(res_deseq_fp, alpha = 0.05)


################################## Normalization Using DESeq2 for Downstream Clustering & Classification

library(vsn)


## We will normalize the data for downstream clustering and classification using DESeq2. Note that applying the DESeq() function normalizes the data internally, and here we are using internal functions to normalize and extract the dataset for downstream use.
## Since our data is in the form of raw RNA sequencing counts, we need to normalize the data by scaling the raw counts to account for differences in counts that is not due to gene expression.
## to do this we can use the median of ratios method of normalization that is included in DESeq2. We can use the estimateSizeFactors() function to generate size factors that we will use to divide the gene count by.
## Normally DESeq2 performs this step automatically, but we will use it here to show the data before and after normalization. 
## SOURCE: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html



# Extract normalized expression matrix with VST, which will normalize the data and compute a variance-stabilizing transformation that is similar to putting the data on the log2 scale. 
# We use blind = TRUE to calculate the across-all-samples variability that will be useful for clustering and classification downstream.
vst_fp <- vst(dds_fp, blind = TRUE)
vst_mat_fp <- assay(vst_fp) # normalized expression matrix

# Check that the variance has been stabilized using meanSdPlot
meanSdPlot(vst_mat_fp)

# convert to expression matrix and phenotype matrix
expr_mat_vst <- vst_mat_fp #normalized expression matrix from DESeq2 with variance stabilized for downstream tasks like clustering.
sample_meta_vst <- as.data.frame(colData(vst_fp)) #phenotype annotation

#check that vst has not removed any samples
dim(expr_mat_vst)
#still 16967 genes and 80 samples.

#### Visualize the normalization done internally via vst using estimateSizeSamples(), which is 
#### applied internally using vst()

dds_fp_vis <- estimateSizeFactors(dds_fp)

# Boxplots of raw counts (before normalization)
par(mfrow = c(1,2))
raw_counts <- counts(dds_fp, normalized = FALSE)
boxplot(log2(raw_counts + 1),
        main = "Before Normalization",
        ylab = "log2(count + 1)",
        col = "lightblue",
        las = 2)

# Boxplot of normalize counts (after DESeq2 normalization)

norm_counts <- counts(dds_fp_vis, normalized = TRUE)
boxplot(log2(norm_counts + 1),
        main = "After Normalization with DESeq2",
        ylab = "log2(count+1)",
        col = "lightblue",
        las = 2)
par(mfrow = c(1,1))

############################################ Gene Set Enrichment Analysis

#### load gene set of interest
library(hypeR) 

HALLMARK <- msigdb_gsets(species = "Homo sapiens", db_species = "HS", collection = "H")

length(HALLMARK$genesets) #50 gene sets in the library

#extract Epithelial Mesenchymal Transition gene set.

gs <- HALLMARK$genesets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]

length(gs) # There are 200 genes in the gene set

# convert deseq2 object to a data frame for gene set enrichment analysis

dds_res_df <- as.data.frame(res_deseq_fp)
dds_res_df$gene_id <- rownames(dds_res_df)
dds_res_df$symbol <- fData(fp_filtered)$hgnc_symbol[match(dds_res_df$gene_id, rownames(fData(fp_filtered)))]

#Convert gene symbols to upper case
gs_upper <- toupper(gs)
symbols_upper <- toupper(dds_res_df$symbol)

#Match gene set to DESeq2 result symbols
gene_idx <- match(gs_upper, symbols_upper)
gene_idx <- gene_idx[!is.na(gene_idx)]
length(gene_idx) #197 genes

# rank genes based on test statistic from DESeq2 ("stat" column)
gene_ranks <- rank(-dds_res_df$stat)[gene_idx] #rank G3 at the top - take genes most upregulated in G3 (highest positive statistic) get ranked highest and first

# check for NAs
stopifnot(!any(is.na(gene_ranks)))

##### Perform Enrichment Testing using the ksGenescore() function from the demo. We will test for enrichment of the gene set with respect to the phenotype "g3" vs "g1"

source("ksGenescore.R")

# calculate K-S gene scores and generate plot
ksGenescore(n.x = nrow(dds_res_df), y = gene_ranks, do.plot = TRUE) #upregulated in g3 versus g1- higher expression in high grade cancer condition

#### test against the P53 gene set

#extract P53 gene set

gs2 <- HALLMARK$genesets[["HALLMARK_P53_PATHWAY"]]

length(gs2) # There are 200 genes in the gene set

#Convert gene symbols to upper case
gs2_upper <- toupper(gs2)
symbols_upper <- toupper(dds_res_df$symbol)

#Match gene set to DESeq2 result symbols
gene_idx2 <- match(gs2_upper, symbols_upper)
gene_idx2 <- gene_idx2[!is.na(gene_idx2)]
length(gene_idx2) #192 genes

# rank genes based on test statistic from DESeq2 ("stat" column)
gene_ranks2 <- rank(-dds_res_df$stat)[gene_idx2]

# check for NAs
stopifnot(!any(is.na(gene_ranks2)))

##### Perform Enrichment Testing using the ksGenescore() function from the demo. 


# calculate K-S gene scores and generate plot
ksGenescore(n.x = nrow(dds_res_df), y = gene_ranks2, do.plot = TRUE) #enriched in g1 vs g3 - higher expression in low grade cancer condition
#upside down but significantly enriched - indicating enrichment in G1 rather than G3


############################## Hierarchical Clustering



library(Biobase)
library(ComplexHeatmap)
library(dplyr)

#check if vst filtered out any genes

expr_mat_vst  #normalized expression matrix from DESeq2
sample_meta_vst  #phenotype annotation

# ensure 'condition' is a character or factor vector
sample_meta_vst$condition <- as.character(sample_meta_vst$condition)

# subset the genes to the genes with top variability to reduce computational load in hierarchical clustering
mad_vals_vst <- apply(expr_mat_vst, 1, mad)
top_genes_vst <- names(sort(mad_vals_vst, decreasing = TRUE))[1:1000]
expr_mat_top_vst <- expr_mat_vst[top_genes_vst,]


# create heatmap annotation
hm_ann_fp_vst <- HeatmapAnnotation(
  condition = sample_meta_vst$condition,
  col = list(condition = c("g1" = "lightblue", "g3" = "darkred"))
)

# hierarchical clustering for rows and columns
# using helper functions provided in the demo
source("hcopt.R")

hc_col_fp_vst <- hcopt(dist(t(expr_mat_top_vst)), method = "ward.D") #columns are samples
hc_row_fp_vst <- hcopt(as.dist(1 - cor(t(expr_mat_top_vst))), method = "ward.D") #rows are genes

# Create heatmap
Heatmap(
  expr_mat_top_vst,
  name = "vst normalized expression",
  top_annotation = hm_ann_fp_vst,
  cluster_rows = hc_row_fp_vst,
  cluster_columns = hc_col_fp_vst,
  column_split = 2,
  row_split = 3,
  show_parent_dend_line = TRUE,
  row_title = "",
  show_column_names = FALSE,
  show_row_names = FALSE
)



# top clustering by sample - most g3 in cluster 2,most g1 in cluster 1
# largely overlap with supervised labels but not perfect. unsupervised can try to recover sample groups but not perfect
# left- clustering in genes. block 1 has lots of high expressed vs sample 2. 
# gene group 1 tend to be high in sample group 1 vs 2, but we don't know those genes
# Next step -> you'd look into the genes in cluster 1 and see which ones to see if they are of interest


################################# Classification Analysis



#### Fit a Naive Bayes classifier using the caret package. Will build a predictive model to classify tumor grade (g1 vs g3) based on gene expression data.

# Train & Test

# As discussed in class, we start by performing a 50/50 train/test split, then 
# putting aside the test set, and building classifiers on the training set.

library(caret)


set.seed(123)
split_idx_vst <- caret::createDataPartition(sample_meta_vst$condition, p = 0.5, list = FALSE) ##split all 80 samples randomly in half using carat.

# train and test for expression data
vst_train_expr <- expr_mat_vst[, split_idx_vst]
vst_test_expr <- expr_mat_vst[, -split_idx_vst]

# train and test for phenotype data
vst_train_pheno <- sample_meta_vst[split_idx_vst, ]
vst_test_pheno <- sample_meta_vst[-split_idx_vst,]


dim(vst_train_expr)
dim(vst_test_expr)
dim(vst_train_pheno)
dim(vst_test_pheno)

# use training set to fit model.
# Fitting a Naive Bayes model

# We start by fitting a single model, without any model selection (i.e. without 
# searching for the best combination of parameters).
getModelInfo(model = "naive_bayes")$naive_bayes$parameters

# show the default grid generation for the NB model
# the function body
getModelInfo(model = "naive_bayes")$naive_bayes$grid

# set it to what you want it
nb_grid_vst <- expand.grid(usekernel = FALSE, laplace = 1, adjust = 1)
nb_grid_vst # a single model


# The function trainControl() is used to specify the type of resampling to 
# perform "model selection" (e.g. cross validation).
# Since we are not searching over any combination of model parameters, we will 
# set it to "none".
nb_control_vst <- trainControl(method = "none")


# Finally, we are ready to fit our Naive Bayes model.

nb_train_vst <- caret::train(
  x = t(vst_train_expr), 
  y = vst_train_pheno$condition, 
  method = "naive_bayes", 
  trControl = nb_control_vst, 
  tuneGrid = nb_grid_vst
)


# look at the classification accuracy of our model.
# apply to training model
nb_train_vst_probs <- caret::extractProb(
  list(model = nb_train_vst)
)
head(nb_train_vst_probs) ## we're looking within the training set

# compare true values to predictions
table(obs = nb_train_vst_probs$obs, pred = nb_train_vst_probs$pred) ##19 observed as g1 were also predicted, 1 misclassified as g3. 19 correctly classified as g3, 1 misclassified as g1 .

# using caret function
nb_train_vst_confusion <- caret::confusionMatrix(
  data = nb_train_vst_probs$pred, 
  reference = nb_train_vst_probs$obs, 
  positive = "g1")

nb_train_vst_confusion
## accuracy is sum of those predicted correctly.

# Next: evaluate how well we actually do on the test set (instead of on the 
# training set)
nb_test_vst <- predict(nb_train_vst, t(vst_test_expr))

nb_test_vst_probs <- caret::extractProb(
  list(model = nb_train_vst), 
  testX = t(vst_test_expr), 
  testY = vst_test_pheno$condition) |> 
  dplyr::filter(dataType == "Test")

nb_test_vst_confusion <- caret::confusionMatrix(
  data = nb_test_vst_probs$pred, 
  reference = nb_test_vst_probs$obs, 
  positive = "g1")

nb_test_vst_confusion


# When estimating the accuracy on the independent test set rather 
# than on the training set, we see reduced accuracy from 95 to 72.5. 
# Accuracy lower on test because test contained initial samples we did not use to train the model.


