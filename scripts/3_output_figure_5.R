### Figure 5 ###

outdir_files="./output_dir/"
indir_files="./output_dir/"

source("./scripts/utils.R")


#################################


### MAYA on ovary ###

{
    sce<-readRDS("./datasets/Ovary/Zhang_sce.rds")
    meta_ovary<-as.data.frame(colData(sce))
    data_ovary<-counts(sce)
    seurat_ovary<-CreateSeuratObject(data_ovary,meta.data=meta_ovary)
    seurat_ovary<-subset(seurat_ovary, subset= treatment_phase=="treatment-naive")
    meta_ovary<-as.data.frame(seurat_ovary@meta.data)
    data_ovary<-GetAssayData(seurat_ovary,slot="count")
    logcpm_ovary<-logcpmNormalization(data_ovary)
    #saveRDS(logcpm_ovary,file.path(outdir_files,"logcpm_ovary.rds"))
    rm("sce")
    meta_ovary$Cell_type<-meta_ovary$cell_subtype
    meta_ovary$Cell_type[which(meta_ovary$cell_type=="EOC")]<-"EOC"
    meta_ovary$Cell_type[which(meta_ovary$cell_subtype %in% c("CAF-1","CAF-2","CAF-3"))]<-"CAF"
    meta_ovary$Cell_type[which(meta_ovary$cell_subtype %in% c("DC-1","DC-2"))]<-"DC"
    #saveRDS(meta_ovary,file.path(outdir_files,"meta_ovary.rds"))
    #seurat
    seurat_ovary <- NormalizeData(seurat_ovary, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_ovary <- FindVariableFeatures(seurat_ovary, selection.method = "vst", nfeatures = 2000)
    seurat_ovary <- ScaleData(seurat_ovary)
    seurat_ovary <- RunPCA(seurat_ovary, features = VariableFeatures(object = seurat_ovary))
    seurat_ovary <- RunUMAP(seurat_ovary, dims = 1:10)
    
    umap_ovary<-Embeddings(seurat_ovary,reduction = "umap")
    umap_ovary<-list(layout=umap_ovary)
    #saveRDS(umap_ovary,file.path(outdir_files,"umap_ovary_gene_based.rds"))
    
    pca_ovary<-Embeddings(seurat_ovary,reduction = "pca")
    
}


# Hallmark
ovary_hallmark<-MAYA_pathway_analysis(expr_mat = data_ovary,
                                      modules_list = "hallmark",
                                      min_cells_pct = 0.01,
                                      is_logcpm = F,
                                      min_genes = 10,
                                      max_contrib = 0.4,
                                      scale_before_pca = T,
                                      all_PCs_in_range = F)

saveRDS(ovary_hallmark,file=file.path(outdir_files,"PCA_obj_hallmark_ovary.rds"))


### Shannon

ovary_patient<-seurat_ovary$patient_id
ovary_patient[which(seurat_ovary$patient_id=="EOC372")]<-1
ovary_patient[which(seurat_ovary$patient_id=="EOC443")]<-2
ovary_patient[which(seurat_ovary$patient_id=="EOC540")]<-3
ovary_patient[which(seurat_ovary$patient_id=="EOC3")]<-4
ovary_patient[which(seurat_ovary$patient_id=="EOC87")]<-5
ovary_patient[which(seurat_ovary$patient_id=="EOC136")]<-6
ovary_patient[which(seurat_ovary$patient_id=="EOC1005")]<-7
ovary_patient[which(seurat_ovary$patient_id=="EOC733")]<-8
ovary_patient[which(seurat_ovary$patient_id=="EOC153")]<-9
ovary_patient[which(seurat_ovary$patient_id=="EOC349")]<-10
ovary_patient[which(seurat_ovary$patient_id=="EOC227")]<-11

# gene-based 
cluster_result_seurat<-cluster_from_pca(pca_ovary,res=0.00005)
df_seurat<-compute_SDI(annotation=ovary_patient,cluster_annot=cluster_result_seurat)

tmp<-table(cluster_result_seurat[which(meta_ovary$Cell_type=="EOC")])
tumor_clust<-names(tmp)[which(tmp>50)]
df_seurat$is_tumor<-ifelse(df_seurat$cluster %in% tumor_clust,TRUE,FALSE)

# MAYA based
cluster_result_maya<-cluster_from_pca(t(ovary_hallmark$activity_matrix),res=0.00005)
df_maya<-compute_SDI(annotation=ovary_patient,cluster_annot=cluster_result_maya)

tmp<-table(cluster_result_maya[which(meta_ovary$Cell_type=="EOC")])
tumor_clust<-names(tmp)[which(tmp>50)]
df_maya$is_tumor<-ifelse(df_maya$cluster %in% tumor_clust,TRUE,FALSE)

final_df_patient<-data.frame(method=c(rep("Gene-based",length(df_seurat$shannon)),rep("MAYA",length(df_maya$shannon))),
                             shannon=c(df_seurat$shannon,df_maya$shannon),
                             is_tumor=c(df_seurat$is_tumor,df_maya$is_tumor))

write.table(final_df_patient,file=file.path(outdir_files,"ovary_shannon_results.tsv"),sep="\t",col.names=T,row.names=F,quote=F)

df<-data.frame(cells=colnames(seurat_ovary),
               gene_based_clusters=paste0("C",cluster_result_seurat),
               maya_clusters=paste0("C",cluster_result_maya))

write.table(df,file=paste0(outdir_files,"ovary_cells_annot_shannon.tsv"),sep="\t",col.names=T,row.names=F,quote=F)

rm(seurat_ovary,logcpm_ovary,data_ovary)


## KEGG

ovary_kegg<-MAYA_pathway_analysis(expr_mat = data_ovary,
                                  modules_list = "kegg",
                                  min_cells_pct = 0.01,
                                  is_logcpm = F,
                                  min_genes = 10,
                                  max_contrib = 0.4)
saveRDS(ovary_kegg,file=file.path(outdir_files,"ovary_kegg.rds"))


#################################




### perform GSEA ###

Idents(seurat_ovary)<-"Cell_type"
seurat_ovary<-NormalizeData(seurat_ovary,normalization.method = "LogNormalize",scale.factor = 10000)

for(cell_type in c("EOC","CAF","Macrophages")){
    logfc<-FoldChange(seurat_ovary,rownames(meta_ovary[which(meta_ovary$Cell_type==cell_type),]),rownames(meta_ovary[which(meta_ovary$Cell_type!=cell_type),]))
    logfc<-logfc[order(logfc$avg_log2FC,decreasing = T),]
    gene_list_DEGs<-logfc[,"avg_log2FC"]
    names(gene_list_DEGs)<-rownames(logfc)
    gsea_res<-GSEA(gene_list_DEGs,hallmark_modules)
    saveRDS(gsea_res,file.path(outdir_files,paste0("gsea_res_",cell_type,".rds")))
}


cell_type="EOC"
logfc<-FoldChange(seurat_ovary,rownames(meta_ovary[which(meta_ovary$Cell_type==cell_type),]),rownames(meta_ovary[which(meta_ovary$Cell_type!=cell_type),]))
logfc<-logfc[order(logfc$avg_log2FC,decreasing = T),]

saveRDS(logfc,file.path(outdir_files,paste0("gsea_res_",cell_type,"_fold_change.rds")))


#################################




### NMF ###

seurat_ovary<-NormalizeData(seurat_ovary,normalization.method = "LogNormalize",scale.factor = 10000)

seurat<-subset(seurat_ovary, subset= cell_type =="EOC")

expr_tumor<-list()
top_var_genes<-list()
for(patient in unique(seurat$patient_id)){
    tmp<-subset(seurat, subset=patient_id==patient)
    expr_tumor[[patient]]<-GetAssayData(tmp,slot="data")
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 7000)
    top_var_genes[[patient]]<-VariableFeatures(tmp)
}

rm(seurat)


# perform NMF with ranks ranging from 4 to 9       
w_basis_tumor <- list() # nmf gene scores
h_coef_tumor <- list() # nmf cell scores


for(i in names(expr_tumor)) {
    w <- NULL
    h <- NULL
    CP100K_log=expr_tumor[[i]]
    CP100K_log <- CP100K_log[top_var_genes[[i]],]
    CP100K_log <- CP100K_log - rowMeans(CP100K_log)
    CP100K_log[CP100K_log < 0] <- 0
    for(j in 4:9) {
        nmf_programs <- RcppML::nmf(CP100K_log, k=j, seed=1)
        nmf_programs_scores <- list(w_basis=nmf_programs$w, h_coef=t(nmf_programs$h))
        colnames(nmf_programs_scores$w_basis) <- paste0(i, "_", j, ".", 1:j)
        colnames(nmf_programs_scores$h_coef) <- paste0(i, "_", j, ".", 1:j)
        
        w <- cbind(w, nmf_programs_scores$w_basis)
        h <- cbind(h, nmf_programs_scores$h_coef)
    }
    rownames(w) <- rownames(CP100K_log)
    w_basis_tumor[[i]] <- w
    rownames(h) <- colnames(CP100K_log)
    h_coef_tumor[[i]] <- h
}

# save output
saveRDS(w_basis_tumor, file.path(outdir_files,"nmf_programs_genes_tumor.rds"))
saveRDS(h_coef_tumor, file.path(outdir_files,"nmf_programs_cells_tumor.rds"))

w_basis_tumor<-readRDS(file.path(outdir_files,"nmf_programs_genes_tumor.rds"))
h_coef_tumor<-readRDS(file.path(outdir_files,"nmf_programs_cells_tumor.rds"))

nmf_programs_genes_tumor <- w_basis_tumor
nmf_programs_cells_tumor <- h_coef_tumor


nmf_programs_sig_tumor <- lapply(nmf_programs_genes_tumor, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))

# for each cell line, select robust NMF programs (i.e. obseved using different ranks in the same cell line), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other cell lines. 
nmf_filter_tumor <- robust_nmf_programs(nmf_programs_sig_tumor, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
nmf_programs_sig_tumor <- lapply(nmf_programs_sig_tumor, function(x) x[, is.element(colnames(x), nmf_filter_tumor),drop=F])
nmf_programs_sig_tumor <- do.call(cbind, nmf_programs_sig_tumor)

# calculate similarity between programs
nmf_intersect_tumor <- apply(nmf_programs_sig_tumor , 2, function(x) apply(nmf_programs_sig_tumor , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_tumor <- hclust(as.dist(50-nmf_intersect_tumor), method="average") 
nmf_intersect_hc_tumor <- reorder(as.dendrogram(nmf_intersect_hc_tumor), colMeans(nmf_intersect_tumor))
nmf_intersect_tumor <- nmf_intersect_tumor[order.dendrogram(nmf_intersect_hc_tumor), order.dendrogram(nmf_intersect_hc_tumor)]

# plot similarity matrix heatmap     
nmf_intersect_melt_tumor <- reshape2::melt(nmf_intersect_tumor) 

saveRDS(nmf_intersect_melt_tumor,file.path(outdir_files,"NMF_correlation_mat.rds"))



# enrichment
gene_sets_msigdb <- load_MSIGdb(ref="hg38",GeneSetClasses=c("hallmark"))

list_enrichment<-list()
program<-c()
pathway<-c()
qvalue<-c()
for(prog in names(NMF_programs)){
    tmp <- enrichment_markers_provided_markers(ref="hg38",GeneSetClasses=c("hallmark"),top_markers = NMF_programs[[prog]],qval = 0.05)
    list_enrichment[[prog]]<-tmp
    if(nrow(tmp)!=0){
        program<-c(program,rep(prog,nrow(tmp)))
        pathway<-c(pathway,rownames(tmp))
        qvalue<-c(qvalue,tmp$`q-value`)
    }
    
}
df<-data.frame(program,pathway,qvalue)
df$log10qval<-(-log10(df$qvalue))

saveRDS(df,file.path(outdir_files,"NMF_enrichment.rds"))


#################################




### Other tumor datasets ###

# breast 

seurat<-readRDS("./datasets/Breast/breast_seurat.rds")

meta_breast<-seurat@meta.data
meta_breast$sample<-as.factor(meta_breast$sample)

logcpm_breast<-logcpmNormalization(GetAssayData(seurat,slot="count"))

breast_hallmark<-MAYA_pathway_analysis(expr_mat = logcpm_breast,
                                       modules_list = "hallmark",
                                       min_cells_pct = 0.01,
                                       is_logcpm = T,
                                       min_genes = 10,
                                       max_contrib = 0.4,
                                       scale_before_pca = T,
                                       all_PCs_in_range = F)

saveRDS(meta_breast,file=file.path(outdir_files,"meta_breast.rds"))
saveRDS(logcpm_breast,file=file.path(outdir_files,"logcpm_breast.rds"))
saveRDS(breast_hallmark,file=file.path(outdir_files,"breast_hallmark.rds"))

rm(seurat)



# lung

seurat<-readRDS("./datasets/Lung/lung_seurat.rds")

meta_lung<-seurat@meta.data
meta_lung$sample<-as.factor(meta_lung$sample)
logcpm_lung<-logcpmNormalization(GetAssayData(seurat,slot="count"))

lung_hallmark<-MAYA_pathway_analysis(expr_mat = logcpm_lung,
                                     modules_list = "hallmark",
                                     min_cells_pct = 0.01,
                                     is_logcpm = T,
                                     min_genes = 10,
                                     max_contrib = 0.4,
                                     scale_before_pca = T,
                                     all_PCs_in_range = F)

saveRDS(meta_lung,file=file.path(outdir_files,"meta_lung.rds"))
saveRDS(logcpm_lung,file=file.path(outdir_files,"logcpm_lung.rds"))
saveRDS(lung_hallmark,file=file.path(outdir_files,"lung_hallmark.rds"))



#################################








