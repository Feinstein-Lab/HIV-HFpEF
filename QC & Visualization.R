library(Seurat)
library(scDblFinder)
###QUALIT CONTROL ON COMBINED SAMPLES#####
# Add in the Mitochondrial PCT% information
Full_srt$percent.mt <- PercentageFeatureSet(Full_srt, pattern = "^MT-")
# nCount_RNA is the number of UMI counts in a cell
hist(Full_srt$nCount_RNA)
# nFeature_RNA is the number of different genes that had any reads
hist(Full_srt$nFeature_RNA)
# percent.mt is the percent mitochondrial reads
hist(Full_srt$percent.mt)

# Make a violin plot of the QC columns
VlnPlot(Full_srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

# Filter out unwanted cells from full carotid srt
Full_srt <- subset(Full_srt, subset = 
                     nFeature_RNA > 200 & nFeature_RNA < 3000 
                   & nCount_RNA > 200 & nCount_RNA < 10000 
                   & percent.mt < 10) 

VlnPlot(Full_srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

### REMOVE DOUBLETS ###
Full_sce <- as.SingleCellExperiment(Full_srt)
Full_sce <- scDblFinder(Full_sce, verbose = TRUE, nfeatures = 1000)
#Total number of Doublets
table(Full_sce$scDblFinder.class)
#Add back into seurat
Full_srt_clean <- as.Seurat(Full_sce)
# Check the number of cell before subset
dim(Full_srt_clean)
# Remove all doublets found
Full_srt_clean <- subset(Full_srt_clean, subset = scDblFinder.class == "singlet")
# Check the number after subset
dim(Full_srt_clean)

### NORMALIZE & SCALE AGAIN ###
# Log-transform the counts
Full_srt_clean <- NormalizeData(Full_srt_clean)
# Find Variable Features
Full_srt_clean <- FindVariableFeatures(Full_srt_clean)
# Scale the data
Full_srt_clean <- ScaleData(Full_srt_clean)
# Run PCA
Full_srt_clean <- RunPCA(Full_srt_clean)
## Choose the number of principle components to keep
ElbowPlot(Full_srt_clean,  ndims = 50)

### DATA VISUALIZATION: PCA ###
#Plot the PCA by sample
DimPlot(Full_srt_clean, reduction = "pca", group.by = "orig.ident")
#change original sample to either "Control" or "Case"
Idents(Full_srt_clean) <- Full_srt_clean$orig.ident
Full_srt_clean <- RenameIdents(Full_srt_clean, 
                               'Control_0134' = 'control',
                               'Control_1603' = 'control',
                               'Control_6147' = 'control',
                               'Control_1363' = 'control',
                               'Case_0299' = 'case',
                               'Case_0921' = 'case',
                               'Case_2568' = 'case',
                               'Case_2759' = 'case'
)

#Add these idents into new column labelled "type"
Full_srt_clean$condition <- Idents(Full_srt_clean)

# Plot the PCA by condition
DimPlot(Full_srt_clean, reduction = "pca", group.by = "condition")

###Determine the % Variation in each PC###
# Access the PCA results
pca_results <- Full_srt_clean[["pca"]]
# Get the standard deviations of the principal components
std_devs <- pca_results@stdev
# Calculate the variance explained by each PC
variances <- (std_devs^2)
# Calculate the total variance
total_variance <- sum(variances)
# Calculate the percentage of variance explained by each PC
percentage_variance_explained <- (variances / total_variance) * 100
# Display the results
print(percentage_variance_explained)
# Create a scree plot
plot(x = 1:length(percentage_variance_explained), y = percentage_variance_explained,
     type = "b", xlab = "PC", ylab = "Percentage of Variance Explained",
     main = "Scree Plot")

### DATA VISUALIZATION: UMAP ###
#Reset Idents to samples
Idents(Full_srt_clean) <- Full_srt_clean$orig.ident
# Run UMAP and add embeddings
Full_srt_clean <- RunUMAP(Full_srt_clean, dims = 1:50)
# Find nearest neighbors and construct the graph
Full_srt_clean <- FindNeighbors(Full_srt_clean, k.param = 50, dims = 1:50)
# Find the clusters
Full_srt_clean <- FindClusters(Full_srt_clean, resolution = 0.25)
# Plot the UMAP with samples
DimPlot(Full_srt_clean, reduction = "umap", group.by = "orig.ident")
# Plot the UMAP with condition
DimPlot(Full_srt_clean, reduction = "umap", group.by = "condition")
# Plot the UMAP with cluster
DimPlot(Full_srt_clean, reduction = "umap", group.by = "seurat_clusters")

#####INTEGRATE WITH SEURAT#####################
# Split the seurat object and integrate with CCA
Full_srt_clean_List <- SplitObject(Full_srt_clean, split.by = "orig.ident")
# Normalize and identify variable features for each dataset independently
Full_srt_clean_List <- lapply(X = Full_srt_clean_List, SCTransform) 
features <- SelectIntegrationFeatures(object.list = Full_srt_clean_List, 
                                      nfeatures = 3000)
Full_srt_clean_List <- PrepSCTIntegration(object.list = Full_srt_clean_List, 
                                          anchor.features = features)
Full_srt_clean_List.anchors <- FindIntegrationAnchors(object.list = Full_srt_clean_List,
                                                      normalization.method = "SCT",
                                                      anchor.features = features)
full_srt.int <- IntegrateData(anchorset = Full_srt_clean_List.anchors, 
                              normalization.method = "SCT")

### NORMALIZE & SCALE AGAIN ###
# Run PCA
full_srt.int <- RunPCA(full_srt.int)
#Find Neighbours
full_srt.int <- FindNeighbors(full_srt.int, k.param = 50, dims = 1:50)
# Find the clusters
full_srt.int <- FindClusters(full_srt.int, resolution = 0.6)
# Get the UMAP embedding
full_srt.int <- RunUMAP(full_srt.int, dims = 1:50)

DimPlot(full_srt.int, reduction = "umap", group.by = "orig.ident")
DimPlot(full_srt.int, reduction = "umap", group.by = "condition")
DimPlot(full_srt.int, reduction = "umap", group.by = "seurat_clusters",label = TRUE) 
DimPlot(full_srt.int, reduction = "umap", split.by = "condition", group.by = "seurat_clusters", label = TRUE)

# Compare no. Cells and proportion of cells in each Cluster
compTable <- table(full_srt.int$condition, full_srt.int$seurat_clusters)
compTable <- (compTable / rowSums(compTable)) * 100
compTable

# Compare no. Cells and proportion of cells in each Cluster
compTable.origident <- table(full_srt.int$orig.ident, full_srt.int$seurat_clusters)
compTable.origident <- (compTable.origident / rowSums(compTable.origident)) * 100
compTable.origident
