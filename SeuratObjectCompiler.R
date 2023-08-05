library(Seurat)

#upload gene expression data
HFpEF0134.data <- Read10X(data.dir = "sample_feature_bc_matrix_0134")
HFpEF0299.data <- Read10X(data.dir = "sample_feature_bc_matrix_0299")
HFpEF0921.data <- Read10X(data.dir = "sample_feature_bc_matrix_0921")
HFpEF1603.data <- Read10X(data.dir = "sample_feature_bc_matrix_1603")
HFpEF2568.data <- Read10X(data.dir = "sample_feature_bc_matrix_2568")
HFpEF2759.data <- Read10X(data.dir = "sample_feature_bc_matrix_2759")
HFpEF6147.data <- Read10X(data.dir = "sample_feature_bc_matrix_6147")
HFpEF1363.data <- Read10X(data.dir = "sample_feature_bc_matrix_1363")

#create Seurat Object
HFpEF0134 <- CreateSeuratObject(counts = HFpEF0134.data, project = "Control_0134")
HFpEF0299 <- CreateSeuratObject(counts = HFpEF0299.data, project = "Case_0299")
HFpEF0921 <- CreateSeuratObject(counts = HFpEF0921.data, project = "Case_0921")
HFpEF1603 <- CreateSeuratObject(counts = HFpEF1603.data, project = "Control_1603")
HFpEF2568 <- CreateSeuratObject(counts = HFpEF2568.data, project = "Case_2568")
HFpEF2759 <- CreateSeuratObject(counts = HFpEF2759.data, project = "Case_2759")
HFpEF6147 <- CreateSeuratObject(counts = HFpEF6147.data, project = "Control_6147")
HFpEF1363 <- CreateSeuratObject(counts = HFpEF1363.data, project = "Control_1363")

Full_srt <- merge(x = HFpEF0134, y = c(HFpEF0299, HFpEF0921, HFpEF1363, HFpEF1603, HFpEF2568, HFpEF2759, HFpEF6147), 
                  add.cell.ids = c("0134", "0299", "0921", "1363","1603", "2568", "2759", "6147"))

#upload TCR contig data
vdj_0134 <- read.csv("filtered_contig_annotations_0134_new.csv")
vdj_0299 <- read.csv("filtered_contig_annotations_0299_new.csv")
vdj_0921 <- read.csv("filtered_contig_annotations_0921_new.csv")
vdj_1603 <- read.csv("filtered_contig_annotations_1603_new.csv")
vdj_2568 <- read.csv("filtered_contig_annotations_2568_new.csv")
vdj_2759 <- read.csv("filtered_contig_annotations_2759_new.csv")
vdj_6147 <- read.csv("filtered_contig_annotations_6147_new.csv")
vdj_1363 <- read.csv("filtered_contig_annotations_1363_new.csv")

# Function to remove the -1 at the end of each barcode. Subsets so only the first line of each barcode is kept as each entry for given barcode will have same clonotype.
barcoder <- function(df, prefix){
  df$barcode_og <- df$barcode
  df <- df[!duplicated(df$barcode), ]
  df$barcode <- paste0(prefix, df$barcode)
  df <- df[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene", "reads", "reads_normalized", "barcode_og", "cdr3")]
  names(df)[names(df) == "raw_clonotype_id"] <- "clonotype_id"
  return(df)
}

# Apply the function to each dataframe 
vdj_0134 <- barcoder(vdj_0134, prefix = "0134_")
vdj_0299 <- barcoder(vdj_0299, prefix = "0299_")
vdj_0921 <- barcoder(vdj_0921, prefix = "0921_")
vdj_1603 <- barcoder(vdj_1603, prefix = "1603_")
vdj_2568 <- barcoder(vdj_2568, prefix = "2568_")
vdj_2759 <- barcoder(vdj_2759, prefix = "2759_")
vdj_6147 <- barcoder(vdj_6147, prefix = "6147_")
vdj_1363 <- barcoder(vdj_1363, prefix = "1363_")

#upload TCR clonotype data
clonotype_0134 <- read.csv("clonotypes_0134.csv")
clonotype_0299 <- read.csv("clonotypes_0299.csv")
clonotype_0921 <- read.csv("clonotypes_0921.csv")
clonotype_1603 <- read.csv("clonotypes_1603.csv")
clonotype_2568 <- read.csv("clonotypes_2568.csv")
clonotype_2759 <- read.csv("clonotypes_2759.csv")
clonotype_6147 <- read.csv("clonotypes_6147.csv")
clonotype_1363 <- read.csv("clonotypes_1363.csv")

# Slap the clonotype data onto our VDJ table by clonotype_id.
vdj_0134 <- merge(vdj_0134, clonotype_0134[, c("clonotype_id", "cdr3s_aa")])
vdj_0299 <- merge(vdj_0299, clonotype_0299[, c("clonotype_id", "cdr3s_aa")])
vdj_0921 <- merge(vdj_0921, clonotype_0921[, c("clonotype_id", "cdr3s_aa")])
vdj_1603 <- merge(vdj_1603, clonotype_1603[, c("clonotype_id", "cdr3s_aa")])
vdj_2568 <- merge(vdj_2568, clonotype_2568[, c("clonotype_id", "cdr3s_aa")])
vdj_2759 <- merge(vdj_2759, clonotype_2759[, c("clonotype_id", "cdr3s_aa")])
vdj_6147 <- merge(vdj_6147, clonotype_6147[, c("clonotype_id", "cdr3s_aa")])
vdj_1363 <- merge(vdj_1363, clonotype_1363[, c("clonotype_id", "cdr3s_aa")])

#combine all vdj data
all_vdj <- rbind(vdj_0134,vdj_0299,vdj_0921,vdj_1603,vdj_2568,vdj_2759,vdj_6147,vdj_1363)
# Reorder so barcodes are first column and set them as rownames.
all_vdj <- all_vdj[, c(2,1,3:11)]
rownames(all_vdj) <- all_vdj$barcode

# Add to the Seurat object's metadata.
Full_srt  <- AddMetaData(object=Full_srt, metadata=all_vdj)
