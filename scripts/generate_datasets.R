### Script to generate input matrices ###

setwd("./datasets/")

library(Seurat)
library(SingleCellExperiment)



### Generate SCExperiment object for Ovary dataset ###

# download GSE165897_UMIcounts_HGSOC.tsv.gz and GSE165897_cellInfo_HGSOC.tsv.gz at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897 and unzip

# rds<-readRDS(file="./Ovary/scexp.RDS")
# mat<-assay(rds)
mat<-read.table(file = "./Ovary/GSE165897_UMIcounts_HGSOC.tsv",header = T,sep="\t")

meta<-read.table(file = "./Ovary/GSE165897_cellInfo_HGSOC.tsv",header = T,sep="\t")
rownames(meta)<-meta$cell

seurat<-CreateSeuratObject(mat,meta.data = meta)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- subset(seurat, subset = nFeature_RNA > 1000 & percent.mt < 20 & nCount_RNA>0 & nCount_RNA<80000)

sce<-as.SingleCellExperiment(merged)
logcounts(sce)<-NULL
saveRDS(sce,"./Ovary/Zhang_sce.rds")

###



### Generate SCExperiment object for Larynx dataset ###

# download GSM4546857_LSCC01_DBEC_UMI.csv at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4546857
# download GSM4546858_LSCC02_DBEC_UMI.csv at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4546858

tissue1<-read.csv("./Larynx/GSM4546857_LSCC01_DBEC_UMI.csv",header = T,skip = 6,check.names = FALSE)
tissue2<-read.csv("./Larynx/GSM4546858_LSCC02_DBEC_UMI.csv",header = T)

rownames(tissue1)<-tissue1$Cell_Index
rownames(tissue2)<-tissue2$X

tissue1<-tissue1[,2:ncol(tissue1)]
tissue2<-tissue2[,2:ncol(tissue2)]

tissue1<-as.matrix(tissue1)
tissue1<-t(tissue1)
tissue2<-as.matrix(tissue2)

tissue1 <- as(tissue1, "dgCMatrix")
tissue2 <- as(tissue2, "dgCMatrix")

# seurat tissue 1
seurat_1 <- CreateSeuratObject(counts = tissue1, project = "tissue1", min.cells = 1)
seurat_1<-RenameCells(seurat_1, add.cell.id = "tissue1")
seurat_1[["percent.mt"]] <- PercentageFeatureSet(seurat_1, pattern = "^MT-")
seurat_1<-subset(seurat_1, subset = nFeature_RNA > 1000 & nCount_RNA < 80000 & percent.mt < 20)
seurat_1$orig.ident<-"tissue1"
seurat_1$patient_id<-"P1"

# seurat tissue 2
seurat_2 <- CreateSeuratObject(counts = tissue2, project = "tissue2", min.cells = 1)
seurat_2<-RenameCells(seurat_2, add.cell.id = "tissue2")
seurat_2[["percent.mt"]] <- PercentageFeatureSet(seurat_2, pattern = "^MT-")
seurat_2<-subset(seurat_2, subset = nFeature_RNA > 1000 & nCount_RNA < 80000 & percent.mt < 20)
seurat_2$orig.ident<-"tissue2"
seurat_2$patient_id<-"P2"

# merge
merged<-merge(seurat_1, y = seurat_2)

sce<-as.SingleCellExperiment(merged)
logcounts(sce)<-NULL
saveRDS(sce,"./Larynx/Song_sce.rds")

###



### Generate SCExperiment object for Breast dataset ###

# download https://www.weizmann.ac.il/sites/3CA/breast : Qian et al.

# put downloaded files in datasets/Breast directory

mat<-readMM("./datasets/Breast/Exp_data_UMIcounts.mtx")
genes<-read.table("./datasets/Breast/Genes.txt",header = F)
cells_info<-read.table("./datasets/Breast/Cells.csv",header = T,sep=",")
colnames(mat)<-cells_info$cell_name
rownames(mat)<-genes$V1
rownames(cells_info)<-cells_info$cell_name

seurat<-CreateSeuratObject(mat,meta.data = cells_info)
saveRDS(seurat,"./datasets/Breast/breast_seurat.rds")


### Generate SCExperiment object for Lung dataset ###

# download https://www.weizmann.ac.il/sites/3CA/lung : Kim et al. 

# put downloaded files in datasets/Lung directory

mat<-readMM("./datasets/Lung/Exp_data_UMIcounts.mtx")
genes<-read.table("./datasets/Lung/Genes.txt",header = F)
cells_info<-read.table("./datasets/Lung/Cells.csv",header = T,sep=",")
colnames(mat)<-cells_info$cell_name
rownames(mat)<-genes$V1
rownames(cells_info)<-cells_info$cell_name

seurat<-CreateSeuratObject(mat,meta.data = cells_info)
saveRDS(seurat,"./datasets/Lung/lung_seurat.rds")




### The SCE objects for the Kidney and Colon datasets are already provided in the repo but you can also use this code to reproduce them. 



### Generate SCExperiment object for Kidney dataset ###

# Download data included as supplementary files in the paper Young et al. 

data <- readMM("./Kidney/matrix.mtx")
cells <- read.table("./Kidney/barcodes.tsv",sep="\t",header=T)
genes <- read.table("./Kidney/genes.tsv",sep="\t",header=T)

colnames(data)<-cells$DropletID
rownames(data)<-genes$GeneLabel

meta<-read.table("./Kidney/table11_cell_manifest.txt",header=T,sep="\t")
meta<-meta[,c("DropletID","ClusterID","Compartment","nUMI","nGenes","MTfrac","QCpass","Source")]

cluster_info<-read.table("../table2_cluster_info.txt",header=T,sep="\t")
#keep ClusterID, Alias, Category, Cell_type1, Cell_type2, Cell_type3, Genotype
cluster_info<-cluster_info[,c("Cluster_ID","Alias","Category","Cell_type1","Cell_type2","Cell_type3","Genotype")]

master_meta<-merge(meta,cluster_info,by.x="ClusterID",by.y="Cluster_ID")
rownames(master_meta)<-master_meta$DropletID
master_meta<-master_meta[,c("ClusterID","Alias","Category","Cell_type1","Cell_type2")]

#loading annotations to select only protein coding genes
load("./annotation/Gencode_hg38_v25.RData")
#select protein coding genes in GencodeByGene
coding<-GencodeByGene[which(GencodeByGene$Gene_biotype=="protein_coding"),"Gene_name"]#19899 
#get geneLabels for those genes
rows_to_keep<-genes[which(genes$Symbol%in%coding),c("GeneLabel","Symbol")]

#subset matrix and change row names to symbols
data<-data[rows_to_keep$GeneLabel,]
rownames(data)<-rows_to_keep$Symbol

normal_cells<-as.character(rownames(master_meta[which(master_meta$Category %in% c("Normal_mature_kidney","Normal_mature_kidney_immune")),]))
normal_data<-data[,normal_cells]
normal_meta<-master_meta[normal_cells,]

AV2<-rownames(normal_meta[which(normal_meta$Alias=="AV2"),])
G<-rownames(normal_meta[which(normal_meta$Alias=="G"),])
T8<-rownames(normal_meta[which(normal_meta$Alias=="8T"),])
MNP1<-rownames(normal_meta[which(normal_meta$Alias=="MNP1"),])
M<-rownames(normal_meta[which(normal_meta$Alias=="M"),])

final_cells<-c(AV2,G,T8,MNP1,M)

#Seurat object
seurat_final<-CreateSeuratObject(normal_data[,final_cells],meta.data=normal_meta[final_cells,])
#saveRDS(seurat_final,"./Kidney/seurat_5_clusters.rds")

sce<-as.SingleCellExperiment(seurat_final)
logcounts(sce)<-NULL
saveRDS(sce,"./Kidney/sce_kidney.rds")

###


### Generate SCExperiment object for Colon dataset ###

# Download GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt and GSE144735_processed_KUL3_CRC_10X_annotation.txt at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735

data <- read.table("./Colon/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt")
colnames(data)<-data[1,]
data<-data[-1,]
rownames(data)<-data[,1]
data<-data[,-1]

#loading annotations to select only protein coding genes
load("C./annotation/Gencode_hg38_v25.RData")
#select protein coding genes in GencodeByGene
coding<-GencodeByGene[which(GencodeByGene$Gene_biotype=="protein_coding"),"Gene_name"]#19899 
#get geneLabels for those genes
genes<-rownames(data)#33694
rows_to_keep<-intersect(genes,coding)#19697

#subset matrix and change row names to symbols
data<-data[rows_to_keep,]

# Loading metadata
meta <- read.table("./Colon/GSE144735_processed_KUL3_CRC_10X_annotation.txt",header=T,sep="\t")
rownames(meta)<-meta$Index
meta<-meta[,-1]

normal_cells<-rownames(meta[which(meta$Class=="Normal"),])
normal<-data[,normal_cells]

cells_subset<-rownames(meta[which(meta$Class=="Normal" & meta$Cell_subtype %in% c("NK cells","Regulatory T cells","Pericytes","Smooth muscle cells","cDC","Proliferating","Mast cells","Goblet cells","Mature Enterocytes","CD19+CD20+ B")),])

seurat_sub<-CreateSeuratObject(normal[,cells_subset],meta.data=meta[cells_subset,])
sce<-as.SingleCellExperiment(seurat_sub)
logcounts(sce)<-NULL
saveRDS(sce,"./Colon/sce_colon.rds")


###



