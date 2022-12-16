### Figure 4 ###

outdir_files="./output_dir/"
indir_files="./output_dir/"

source("./scripts/utils.R")


#################################

### Identity ###


# Load data

kidney_list<-generate_panglao_list(organs="Kidney")
digestive_list<-generate_panglao_list(organs="GI tract")

logcpm_kidney<-readRDS(file.path(indir_files,"logcpm_kidney.rds"))
data_kidney<-readRDS(file.path(indir_files,"data_kidney.rds"))
seurat_kidney<-readRDS(file.path(indir_files,"seurat_kidney.rds"))
meta_kidney<-readRDS(file.path(indir_files,"meta_kidney.rds"))

logcpm_crc<-readRDS(file.path(indir_files,"logcpm_crc.rds"))
data_crc<-readRDS(file.path(indir_files,"data_crc.rds"))
seurat_crc<-readRDS(file.path(indir_files,"seurat_crc.rds"))
meta_crc<-readRDS(file.path(indir_files,"meta_crc.rds"))


### MAYA

kidney_panglao<-MAYA_predict_cell_types(expr_mat = data_kidney,
                                     modules_list = kidney_list,
                                     min_cells_pct = 0.05,
                                     organs = NULL,
                                     is_logcpm = FALSE,
                                     nCores = 1,
                                     thr = 0)
saveRDS(kidney_panglao,file.path(outdir_files,"kidney_panglao.rds"))
df<-data.frame(cells=colnames(data_kidney),cell_type=kidney_panglao$cell_annotation)
write.table(df,file.path(outdir_files,"kidney_MAYA_cell_annot.tsv"),col.names = T,row.names = F,quote=F,sep="\t")


crc_panglao<-MAYA_predict_cell_types(expr_mat = data_crc,
                                  modules_list = digestive_list,
                                  min_cells_pct = 0.01,
                                  organs = NULL,
                                  is_logcpm = FALSE,
                                  nCores = 1,
                                  thr=0)
saveRDS(crc_panglao,file.path(outdir_files,"crc_panglao.rds"))
df<-data.frame(cells=colnames(data_crc),cell_type=crc_panglao$cell_annotation)
write.table(df,file.path(outdir_files,"CRC_MAYA_cell_annot.tsv"),col.names = T,row.names = F,quote=F,sep="\t")



### CellID

# kidney

kidney_signature <- lapply(kidney_list, function(x) intersect(x,rownames(seurat_kidney)))
Idents(seurat_kidney)<-"Cell_type"
seurat_kidney <- RunMCA(seurat_kidney, nmcs = 50)
{
    set.seed(1)
    cite_HGT <-
        RunCellHGT(
            X = seurat_kidney,
            reduction = "mca",
            pathways = kidney_signature,
            log.trans = T, minSize = 5
        )
    write.table(as.matrix(cite_HGT),file=paste0(outdir_files,"kidney","_CellID.tsv"),col.names = T,row.names = T,quote=F)
    pred <- rownames(cite_HGT)[apply(cite_HGT, 2, which.max)]
    pred_CellID <-ifelse(2 < apply(cite_HGT, 2, max), pred, "unassigned")
    df<-data.frame(cells=colnames(seurat_kidney),cell_type=pred_CellID)
    write.table(df,file=paste0(outdir_files,"kidney","_CellID_cell_annot.tsv"),col.names = T,row.names = F,quote=F,sep="\t")
    
}

# CRC

digestive_signature <- lapply(digestive_list, function(x) intersect(x,rownames(seurat_crc)))
Idents(seurat_crc)<-"Cell_type"
seurat_crc <- RunMCA(seurat_crc, nmcs = 50)
{
    set.seed(1)
    cite_HGT <-
        RunCellHGT(
            X = seurat_crc,
            reduction = "mca",
            pathways = digestive_signature,
            log.trans = T, minSize = 5
        )
    write.table(as.matrix(cite_HGT),file=paste0(outdir_files,"CRC","_CellID.tsv"),col.names = T,row.names = T,quote=F)
    pred <- rownames(cite_HGT)[apply(cite_HGT, 2, which.max)]
    pred_CellID <-ifelse(2 < apply(cite_HGT, 2, max), pred, "unassigned")
    df<-data.frame(cells=colnames(seurat_crc),cell_type=pred_CellID)
    write.table(df,file=paste0(outdir_files,"CRC","_CellID_cell_annot.tsv"),col.names = T,row.names = F,quote=F,sep="\t")
    
}

rm(seurat_kidney,seurat_crc)


### SCINA

# Kidney
kidney_signature <- lapply(kidney_list, function(x) intersect(x,rownames(logcpm_kidney)))

results = SCINA(exp = logcpm_kidney, signatures = kidney_signature, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE, log_file=paste0(outdir_files,'SCINA_kidney.log'))
df<-data.frame(cells=colnames(logcpm_kidney),cell_type=results$cell_labels)
write.table(df,file=paste0(outdir_files,"kidney","_SCINA_cell_annot.tsv"),col.names = T,row.names = F,quote=F,sep="\t")


# CRC
digestive_signature <- lapply(digestive_list, function(x) intersect(x,rownames(logcpm_kidney)))

results = SCINA(exp = logcpm_crc, signatures = digestive_signature, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE, log_file=paste0(outdir_files,'SCINA_crc.log'))
df<-data.frame(cells=colnames(logcpm_crc),cell_type=results$cell_labels)
write.table(df,file=paste0(outdir_files,"CRC","_SCINA_cell_annot.tsv"),col.names = T,row.names = F,quote=F,sep="\t")



### Compute F1-scores

# kidney

map_cite_kidney <- list(MNP1=c("Monocytes"),
                        Endothelium=c("Endothelial cells"),
                        Mesangial_cells=c("Mesangial cells","Smooth muscle cells"),
                        Podocytes=c("Podocytes"),
                        TCD8=c("T cells","T memory cells", "T cytotoxic cells") )


meta_kidney$MAYA<-read.table(file.path(outdir_files,"kidney_MAYA_cell_annot.tsv"),header = T,check.names = F,sep="\t")$cell_type
meta_kidney$CellID<-read.table(file.path(outdir_files,"kidney_CellID_cell_annot.tsv"),header = T,check.names = F,sep="\t")$cell_type
meta_kidney$SCINA<-read.table(file.path(outdir_files,"kidney_SCINA_cell_annot.tsv"),header = T,check.names = F,sep="\t")$cell_type

F1_list<-c()
for(method in c("CellID","SCINA","MAYA")){
    F1_list[[method]]<-compute_F1_score(predictions=method,reference="Cell_type",meta = meta_kidney,mapping = map_cite_kidney,data =method)
}

df<-rbind(F1_list[[1]],F1_list[[2]],F1_list[[3]])

write.table(df,file=file.path(outdir_files,"kidney_F1_scores.tsv"),sep="\t",quote=F,col.names=T,row.names = F)



# CRC

map_cite_CRC<- list(`Mature Enterocytes`=c("Enterocytes"),
                    `Goblet cells`=c("Goblet cells"),
                    Pericytes=c("Pericytes"),
                    `Smooth muscle cells`=c("Smooth muscle cells"),
                    cDC=c("Dendritic cells"),
                    Proliferating=c("Monocytes","Macrophages"),
                    `NK cells`=c("NK cells","Natural killer T cells"),
                    `Regulatory T cells`=c("T regulatory cells","T cells","T memory cells"),
                    `CD19+CD20+ B`=c("B cells","B cells naive","B cells memory"),
                    `Mast cells`=c("Mast cells"))

meta_crc$MAYA<-read.table(file.path(outdir_files,"CRC_MAYA_cell_annot.tsv"),header = T,check.names = F,sep="\t")$cell_type
meta_crc$CellID<-read.table(file.path(outdir_files,"CRC_CellID_cell_annot.tsv"),header = T,check.names = F,sep="\t")$cell_type
meta_crc$SCINA<-read.table(file.path(outdir_files,"CRC_SCINA_cell_annot.tsv"),header = T,check.names = F,sep="\t")$cell_type

F1_list<-c()
for(method in c("CellID","SCINA","MAYA")){
    F1_list[[method]]<-compute_F1_score(predictions=method,reference="Cell_subtype",meta = meta_crc,mapping = map_cite_CRC,data =method)
}

df<-rbind(F1_list[[1]],F1_list[[2]],F1_list[[3]])

write.table(df,file=file.path(outdir_files,"CRC_F1_scores.tsv"),sep="\t",quote=F,col.names=T,row.names = F)


#################################



### Ovary confusion matrix ###

{
    sce<-readRDS("./datasets/Ovary/Zhang_sce.rds")
    meta_ovary<-as.data.frame(colData(sce))
    data_ovary<-counts(sce)
    seurat_ovary<-CreateSeuratObject(data_ovary,meta.data=meta_ovary)
    seurat_ovary<-subset(seurat_ovary, subset= treatment_phase=="treatment-naive")
    meta_ovary<-as.data.frame(seurat_ovary@meta.data)
    data_ovary<-GetAssayData(seurat_ovary,slot="count")
    logcpm_ovary<-logcpmNormalization(data_ovary)
    saveRDS(logcpm_ovary,file.path(outdir_files,"logcpm_ovary.rds"))
    rm("sce")
    meta_ovary$Cell_type<-meta_ovary$cell_subtype
    meta_ovary$Cell_type[which(meta_ovary$cell_type=="EOC")]<-"EOC"
    meta_ovary$Cell_type[which(meta_ovary$cell_subtype %in% c("CAF-1","CAF-2","CAF-3"))]<-"CAF"
    meta_ovary$Cell_type[which(meta_ovary$cell_subtype %in% c("DC-1","DC-2"))]<-"DC"
    saveRDS(meta_ovary,file.path(outdir_files,"meta_ovary.rds"))
    #seurat
    seurat_ovary <- NormalizeData(seurat_ovary, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_ovary <- FindVariableFeatures(seurat_ovary, selection.method = "vst", nfeatures = 2000)
    seurat_ovary <- ScaleData(seurat_ovary)
    seurat_ovary <- RunPCA(seurat_ovary, features = VariableFeatures(object = seurat_ovary))
    seurat_ovary <- RunUMAP(seurat_ovary, dims = 1:10)
    
    umap_ovary<-Embeddings(seurat_ovary,reduction = "umap")
    umap_ovary<-list(layout=umap_ovary)
    saveRDS(umap_ovary,file.path(outdir_files,"umap_ovary_gene_based.rds"))
    
    pca_ovary<-Embeddings(seurat_ovary,reduction = "pca")
    
}

# Panglao
ovary_panglao<-MAYA_predict_cell_types(expr_mat = data_ovary,
                                    modules_list = NULL,
                                    min_cells_pct = 0.001,
                                    organs = "Ovary",
                                    is_logcpm = FALSE,
                                    nCores = 1)

saveRDS(ovary_panglao,file=file.path(outdir_files,"ovary_panglao.rds"))


#################################



### Batch effect ###

# Larynx 

# Loading dataset
{
    sce<-readRDS("./datasets/Larynx/Song_sce.rds")
    data_larynx<-counts(sce)
    meta_larynx<-as.data.frame(colData(sce))
    rm("sce")
    
    # Create Seurat object and process with usual parameters
    seurat_larynx<-CreateSeuratObject(data_larynx,meta.data = meta_larynx)
    seurat_larynx <- NormalizeData(seurat_larynx, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_larynx <- FindVariableFeatures(seurat_larynx, selection.method = "vst", nfeatures = 2000)
    seurat_larynx <- ScaleData(seurat_larynx)
    seurat_larynx <- RunPCA(seurat_larynx, features = VariableFeatures(object = seurat_larynx))
    seurat_larynx <- RunUMAP(seurat_larynx, dims = 1:10)
    
    umap_larynx<-Embeddings(seurat_larynx,reduction = "umap")
    umap_larynx<-list(layout=umap_larynx)
    saveRDS(umap_larynx,file.path(indir_files,"UMAP_larynx_gene_based.rds"))
    
    logcpm_larynx<-logcpmNormalization(data_larynx)
    
    pca_larynx<-Embeddings(seurat_larynx,reduction = "pca")
}

# Annotate each dataset with cell types using markers given by author
{
    P1<-subset(seurat_larynx,subset= patient_id == "P1")
    P1 <- NormalizeData(P1, normalization.method = "LogNormalize", scale.factor = 10000)
    
    #findvariablefeatures
    P1 <- FindVariableFeatures(P1, selection.method = "vst", nfeatures = 2000)
    
    #scaling data
    P1 <- ScaleData(P1) # remove , features = all.genes
    P1 <- RunPCA(P1, features = VariableFeatures(object = P1))
    
    P1 <- FindNeighbors(P1, dims = 1:10)
    P1 <- FindClusters(P1, resolution = 0.1)
    P1 <- FindClusters(P1, resolution = 1)
    
    P1 <- RunUMAP(P1, dims = 1:10)
    
    P1_type<-rep("unknown",length(P1$RNA_snn_res.0.1))
    P1_type[which(P1$RNA_snn_res.0.1==0 | P1$RNA_snn_res.0.1==2)]<-"Tumor_cells"
    P1_type[which(P1$RNA_snn_res.0.1==1)]<-"T_cells"
    P1_type[which(P1$RNA_snn_res.0.1==3)]<-"Macrophages"
    P1_type[which(P1$RNA_snn_res.0.1==4)]<-"Fibroblasts"
    P1_type[which(P1$RNA_snn_res.0.1==5)]<-"Epithelial_cells"
    P1_type[which(P1$RNA_snn_res.0.1==6)]<-"Endothelial_cells"
    
    P2<-subset(seurat_larynx,subset= patient_id == "P2")
    P2 <- NormalizeData(P2, normalization.method = "LogNormalize", scale.factor = 10000)
    
    #findvariablefeatures
    P2 <- FindVariableFeatures(P2, selection.method = "vst", nfeatures = 2000)
    
    #scaling data
    P2 <- ScaleData(P2) # remove , features = all.genes
    P2 <- RunPCA(P2, features = VariableFeatures(object = P1))
    
    P2 <- FindNeighbors(P2, dims = 1:10)
    P2 <- FindClusters(P2, resolution = 0.1)
    P2 <- FindClusters(P2, resolution = 1)
    
    P2 <- RunUMAP(P2, dims = 1:10)
    
    P2_type<-rep("unknown",length(P2$RNA_snn_res.0.1))
    P2_type[which(P2$RNA_snn_res.0.1==0 | P2$RNA_snn_res.0.1==3)]<-"Tumor_cells"
    P2_type[which(P2$RNA_snn_res.0.1==1)]<-"T_cells"
    P2_type[which(P2$RNA_snn_res.0.1==2)]<-"Macrophages"
    P2_type[which(P2$RNA_snn_res.0.1==5)]<-"Fibroblasts"
    P2_type[which(P2$RNA_snn_res.0.1==4)]<-"Endothelial_cells"
    
    meta_larynx$Cell_type<-c(P1_type,P2_type)
}
saveRDS(meta_larynx,file.path(outdir_files,"meta_larynx.rds"))

# Run MAYA

test_larynx<-MAYA_predict_cell_types(expr_mat = data_larynx,
                                   modules_list = NULL,
                                   min_cells_pct = 0.01,
                                   organs = NULL,
                                   is_logcpm = FALSE,
                                   nCores = 1)

# save annotation and MAYA object
write.table(data.frame(cells=colnames(seurat_larynx),annotation=test_larynx$cell_annotation),file=file.path(outdir_files,"larynx_maya_annotation.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
saveRDS(test_larynx, file=file.path(outdir_files,"larynx_panglao.rds"))



### Compute Shannon diversity index ###

larynx_patient<-seurat_larynx$patient_id
larynx_patient[which(seurat_larynx$patient_id=="P1")]<-1
larynx_patient[which(seurat_larynx$patient_id=="P2")]<-2

# Seurat
cluster_result_seurat<-cluster_from_pca(pca_larynx,res=0.001)
df_seurat<-compute_SDI(annotation=larynx_patient,cluster_annot=cluster_result_seurat)

# Harmony

{
    seurat_larynx <- seurat_larynx %>%
        RunHarmony("patient_id", plot_convergence = TRUE)
    seurat_larynx<- seurat_larynx %>%
        RunUMAP(reduction = "harmony", dims = 1:10)
    
    pca_harmony<-Embeddings(seurat_larynx,reduction = "harmony")
    
    umap_harmony<-Embeddings(seurat_larynx,reduction = "umap")
    umap_harmony<-list(layout=umap_harmony)
}
saveRDS(umap_harmony, file=file.path(outdir_files,"umap_larynx_harmony.rds"))

cluster_result_harmony<-cluster_from_pca(pca_harmony,res=0.001)
df_harmony<-compute_SDI(annotation=larynx_patient,cluster_annot=cluster_result_harmony)


# Seurat CCA

{
    integrated<-run_seuratCCA(seurat_larynx,annot_name="patient_id")
    # Check out integration
    UMAPPlot(integrated,group.by="patient_id")
    
    umap_anchors<-Embeddings(integrated,reduction = "umap")
    umap_anchors<-list(layout=umap_anchors)
    
    pca_anchors<-Embeddings(integrated,reduction = "pca")
}

saveRDS(umap_anchors,file=file.path(outdir_files,"larynx_anchors_UMAP_embeddings_seuratCCA.rds"))
saveRDS(pca_anchors,file=file.path(outdir_files,"larynx_anchors_PCA_embeddings_seuratCCA.rds"))

cluster_result_anchors<-cluster_from_pca(pca_anchors,res=0.0005)
df_anchors<-compute_SDI(annotation=larynx_patient,cluster_annot=cluster_result_anchors)


# MAYA

cluster_result_maya<-cluster_from_pca(t(test_larynx$activity_matrix),res=0.001)
df_maya<-compute_SDI(annotation=larynx_patient,cluster_annot=cluster_result_maya)



## Dataframe for boxplot

final_df_patient<-data.frame(method=c(rep("Gene-based",length(df_seurat$shannon)),rep("Harmony",length(df_harmony$shannon)),rep("SeuratCCA",length(df_anchors$shannon)),rep("MAYA",length(df_maya$shannon))),
                             shannon=c(df_seurat$shannon,df_harmony$shannon,df_anchors$shannon,df_maya$shannon))

write.table(final_df_patient,file=file.path(outdir_files,"larynx_shannon_results.tsv"),sep="\t",col.names=T,row.names=F,quote=F)

df<-data.frame(cells=colnames(seurat_larynx),
               gene_based_clusters=paste0("C",cluster_result_seurat),
               harmony_clusters=paste0("C",cluster_result_harmony),
               seuratcca_clusters=paste0("C",cluster_result_anchors),
               maya_clusters=paste0("C",cluster_result_maya))

write.table(df,file=file.path(outdir_files,"larynx_cells_annot_clusters_all_methods.tsv"),sep="\t",col.names=T,row.names=F,quote=F)




# Pancreas 


# Loading dataset


{
    data("panc8")
    seurat_panc8 <- panc8
    rm(panc8)
    
    seurat_panc8 <- NormalizeData(seurat_panc8, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_panc8 <- FindVariableFeatures(seurat_panc8, selection.method = "vst", nfeatures = 2000)
    seurat_panc8 <- ScaleData(seurat_panc8)
    seurat_panc8 <- RunPCA(seurat_panc8, features = VariableFeatures(object = seurat_panc8))
    seurat_panc8 <- RunUMAP(seurat_panc8, dims = 1:10)
    
    umap_panc8<-Embeddings(seurat_panc8,reduction = "umap")
    umap_panc8<-list(layout=umap_panc8)
    saveRDS(umap_panc8,file.path(outdir_files,"umap_pancreas_gene_based.rds"))
    
    data_panc8<-GetAssayData(seurat_panc8,slot="count")
    logcpm_panc8<-logcpmNormalization(data_panc8)
    
    pca_panc8<-Embeddings(seurat_panc8,reduction = "pca")
    
    meta_panc8<-seurat_panc8@meta.data
    saveRDS(meta_panc8,file.path(outdir_files,"meta_panc8.rds"))
}

# Run MAYA

seurat_panc8 <- NormalizeData(seurat_panc8, normalization.method = "RC", scale.factor = 10000)

test_panc8<-MAYA_predict_cell_types(expr_mat = GetAssayData(seurat_panc8,slot = "data"),
                                    modules_list = NULL,
                                    min_cells_pct = 0.0001,
                                    organs = "Pancreas",
                                    is_logcpm = TRUE,
                                    nCores = 1)

# save annotation and MAYA object
write.table(data.frame(cells=colnames(seurat_panc8),annotation=test_panc8$cell_annotation),file=paste0(outdir_files,"pancreas_maya_annotation.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
saveRDS(test_panc8, file=file.path(outdir_files,"pancreas_panglao.rds"))


### Shannon diversity index

panc8_Technology<-seurat_panc8$tech
panc8_Technology[which(seurat_panc8$tech=="celseq")]<-1
panc8_Technology[which(seurat_panc8$tech=="celseq2")]<-2
panc8_Technology[which(seurat_panc8$tech=="fluidigmc1")]<-3
panc8_Technology[which(seurat_panc8$tech=="indrop")]<-4
panc8_Technology[which(seurat_panc8$tech=="smartseq2")]<-5

# Seurat
cluster_result_seurat<-cluster_from_pca(pca_panc8,res=0.0001)
df_seurat<-compute_SDI(annotation=panc8_Technology,cluster_annot=cluster_result_seurat)

# Harmony

{
    seurat_panc8 <- seurat_panc8 %>%
        RunHarmony("tech", plot_convergence = TRUE)
    seurat_panc8<- seurat_panc8 %>%
        RunUMAP(reduction = "harmony", dims = 1:10)
    
    pca_harmony<-Embeddings(seurat_panc8,reduction = "harmony")
    
    umap_harmony<-Embeddings(seurat_panc8,reduction = "umap")
    umap_harmony<-list(layout=umap_harmony)
}
saveRDS(pca_harmony, file=file.path(outdir_files,"pca_panc8_harmony.rds"))
saveRDS(umap_harmony, file=file.path(outdir_files,"umap_panc8_harmony.rds"))

cluster_result_harmony<-cluster_from_pca(pca_harmony,res=0.00005)
df_harmony<-compute_SDI(annotation=panc8_Technology,cluster_annot=cluster_result_harmony)


# Seurat CCA

{
    integrated<-run_seuratCCA(seurat_panc8,annot_name="tech")
    # Check out integration
    UMAPPlot(integrated,group.by="tech")
    
    umap_anchors<-Embeddings(integrated,reduction = "umap")
    umap_anchors<-list(layout=umap_anchors)
    
    pca_anchors<-Embeddings(integrated,reduction = "pca")
}

saveRDS(umap_anchors,file=paste0(outdir_files,"panc8_anchors_UMAP_embeddings_seuratCCA.rds"))
saveRDS(pca_anchors,file=paste0(outdir_files,"panc8_anchors_PCA_embeddings_seuratCCA.rds"))

cluster_result_anchors<-cluster_from_pca(pca_anchors,res=0.0005)
df_anchors<-compute_SDI(annotation=panc8_Technology,cluster_annot=cluster_result_anchors)


# MAYA
cluster_result_maya<-cluster_from_pca(t(test_panc8$activity_matrix),res=0.0002)
df_maya<-compute_SDI(annotation=panc8_Technology,cluster_annot=cluster_result_maya)

# final dataframe
final_df_Technology<-data.frame(method=c(rep("Gene-based",length(df_seurat$shannon)),rep("Harmony",length(df_harmony$shannon)),rep("SeuratCCA",length(df_anchors$shannon)),rep("MAYA",length(df_maya$shannon))),
                                shannon=c(df_seurat$shannon,df_harmony$shannon,df_anchors$shannon,df_maya$shannon))

df<-data.frame(cells=colnames(seurat_panc8),
               gene_based_clusters=paste0("C",cluster_result_seurat),
               harmony_clusters=paste0("C",cluster_result_harmony),
               seuratcca_clusters=paste0("C",cluster_result_anchors),
               maya_clusters=paste0("C",cluster_result_maya))

write.table(df,file=paste0(outdir_files,"panc8_cells_annot_clusters_all_methods.tsv"),sep="\t",col.names=T,row.names=F,quote=F)

write.table(final_df_Technology,file=paste0(outdir_files,"panc8_shannon_results.tsv"),sep="\t",col.names=T,row.names=F,quote=F)



#################################








