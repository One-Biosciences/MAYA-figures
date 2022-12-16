### Figure 1 & 2 & 3 ###

outdir_files="./output_dir/"
dir.create(outdir_files)

### Generate output ###

source("./scripts/utils.R")

hallmark_modules=read_gmt("./datasets/databases/h.all.v7.4.symbols.gmt")
kegg_modules=read_gmt("./datasets/databases/c2.cp.kegg.v7.4.symbols.gmt")
reactome_modules=read_gmt("./datasets/databases/c2.cp.reactome.v7.5.1.symbols.gmt")
reactome_modules=reactome_modules[sapply(reactome_modules,function(x) length(x)>100 & length(x)<300)]


## Running MAYA with HALLMARK list Kidney dataset ##
{
    # Load SingleCellExperiment
    sce_5_clusters<-readRDS("./datasets/Kidney/sce_kidney.rds")
    data_kidney<-counts(sce_5_clusters)
    saveRDS(data_kidney,file.path(outdir_files,"data_kidney.rds"))
    # #meta_data
    meta_kidney<-as.data.frame(colData(sce_5_clusters))
    meta_kidney$Cell_type<-meta_kidney$Cell_type1
    meta_kidney$Cell_type[which(meta_kidney$Cell_type1=="Nephron_epithelium")]<-"Podocytes"
    meta_kidney$Cell_type[which(meta_kidney$Cell_type1=="Nephron_others")]<-"Mesangial_cells"
    meta_kidney$Cell_type[which(meta_kidney$Cell_type1=="CD8 T cell")]<-"TCD8"
    saveRDS(meta_kidney,file.path(outdir_files,"meta_kidney.rds"))
    #seurat basic processing
    seurat_kidney<-CreateSeuratObject(data_kidney,meta.data = meta_kidney)
    seurat_kidney <- NormalizeData(seurat_kidney, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_kidney <- FindVariableFeatures(seurat_kidney, selection.method = "vst", nfeatures = 2000)
    seurat_kidney <- ScaleData(seurat_kidney)
    seurat_kidney <- RunPCA(seurat_kidney, features = VariableFeatures(object = seurat_kidney))
    seurat_kidney <- RunUMAP(seurat_kidney, dims = 1:10)

    umap_kidney<-Embeddings(seurat_kidney,reduction = "umap")
    umap_kidney<-list(layout=umap_kidney)
    saveRDS(umap_kidney,file.path(outdir_files,"umap_kidney_gene_based.rds"))

    logcpm_kidney<-logcpmNormalization(data_kidney)
    saveRDS(logcpm_kidney,file.path(outdir_files,"logcpm_kidney.rds"))
    
    saveRDS(seurat_kidney,file.path(outdir_files,"seurat_kidney.rds"))
    rm(sce_5_clusters,seurat_kidney)
}

kidney_hallmark<-MAYA_pathway_analysis(expr_mat = data_kidney,
                                       modules_list = "hallmark",
                                       min_cells_pct = 0.05,
                                       is_logcpm = F,
                                       min_genes = 10,
                                       max_contrib = 0.4,
                                       scale_before_pca = T,
                                       all_PCs_in_range = F)

saveRDS(kidney_hallmark,file=paste0(outdir_files,"PCA_obj_hallmark_kidney.rds"))


## Running MAYA with KEGG list Colon dataset ##

{
    sce_crc<-readRDS("./datasets/Colon/sce_colon.rds")
    data_crc<-counts(sce_crc)
    saveRDS(data_crc,file.path(outdir_files,"data_crc.rds"))
    meta_crc<-as.data.frame(colData(sce_crc))
    saveRDS(meta_crc,file.path(outdir_files,"meta_crc.rds"))
    #seurat
    seurat_crc<-CreateSeuratObject(data_crc,meta.data = meta_crc)
    seurat_crc <- NormalizeData(seurat_crc, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_crc <- FindVariableFeatures(seurat_crc, selection.method = "vst", nfeatures = 2000)
    seurat_crc <- ScaleData(seurat_crc)
    seurat_crc <- RunPCA(seurat_crc, features = VariableFeatures(object = seurat_crc))
    seurat_crc <- RunUMAP(seurat_crc, dims = 1:10)

    umap_crc<-Embeddings(seurat_crc,reduction = "umap")
    umap_crc<-list(layout=umap_crc)
    saveRDS(umap_crc,file.path(outdir_files,"umap_crc_gene_based.rds"))

    logcpm_crc<-logcpmNormalization(data_crc)
    saveRDS(logcpm_crc,file.path(outdir_files,"logcpm_crc.rds"))
    
    saveRDS(seurat_crc,file.path(outdir_files,"seurat_crc.rds"))
    rm(sce_crc,seurat_crc)
}


crc_kegg<-MAYA_pathway_analysis(expr_mat = data_crc,
                                modules_list = "kegg",
                                min_cells_pct = 0.05,
                                is_logcpm = F,
                                min_genes = 10,
                                max_contrib = 0.4,
                                scale_before_pca = T,
                                all_PCs_in_range = F)

saveRDS(crc_kegg,file=paste0(outdir_files,"PCA_obj_KEGG_crc.rds"))

## Running MAYA with REACTOME list Colon dataset ##

crc_reactome<-MAYA_pathway_analysis(expr_mat = data_crc,
                                    modules_list = reactome_modules,
                                    min_cells_pct = 0.05,
                                    is_logcpm = F,
                                    min_genes = 10,
                                    max_contrib = 0.4,
                                    scale_before_pca = T,
                                    all_PCs_in_range = F)

saveRDS(crc_reactome,file=paste0(outdir_files,"PCA_obj_reactome_crc.rds"))


### Running MAYA with no filters on informativity

# hallmark
no_filt_MAYA_kidney_hallmark<-run_activity_analysis_no_filter_info(expr_mat = logcpm_kidney,
                                                                   modules_list = hallmark_modules,
                                                                   norm = F,
                                                                   nb_comp_max = 5)

# kegg
no_filt_colon_kegg<-run_activity_analysis_no_filter_info(expr_mat = logcpm_crc,
                                                         modules_list = kegg_modules,
                                                         norm = F,
                                                         nb_comp_max = 5)

# reactome
no_filt_colon_reactome<-run_activity_analysis_no_filter_info(expr_mat = logcpm_crc,
                                                             modules_list = reactome_modules,
                                                             norm = F,
                                                             nb_comp_max = 5)


saveRDS(no_filt_MAYA_kidney_hallmark,file=paste0(outdir_files,"PCA_obj_no_info_hallmark_kidney.rds"))
saveRDS(no_filt_colon_kegg,file=paste0(outdir_files,"PCA_obj_no_info_KEGG_crc.rds"))
saveRDS(no_filt_colon_reactome,file=paste0(outdir_files,"PCA_obj_no_info_reactome_crc.rds"))



## Running AUCell and Pagoda2 for comparison ## 

# pagoda2

hallmark.env <- list2env(hallmark_modules)
kegg.env <- list2env(kegg_modules)
reactome.env <- list2env(reactome_modules)

matrices_pagoda2<-list()
names<-c("kidney_hallmark","kidney_kegg","kidney_reactome","crc_hallmark","crc_kegg","crc_reactome")
i=1
for(counts in c(data_kidney,data_crc)){
    for(env in c(hallmark.env,kegg.env,reactome.env)){
        
        counts <- counts[rowSums(counts)>=10, ]
        rownames(counts) <- make.unique(rownames(counts))
        r <- Pagoda2$new(counts, log.scale=TRUE, n.cores=1)
        r$adjustVariance(plot=TRUE, gam.k=10)
        
        r$testPathwayOverdispersion(env, verbose=TRUE, correlation.distance.threshold=0.95, recalculate.pca=FALSE, top.aspects=15)
        
        # restricted matrix
        valid_pathways=r$misc$pathwayODInfo$name[r$misc$pathwayODInfo$valid]
        
        mat<-matrix(data = NA,nrow=length(valid_pathways),ncol=ncol(counts))
        #colnames(mat)<-colnames(r$misc[['pwpca']][[1]]$xp$scores)
        colnames(mat)<-colnames(counts)
        rownames(mat)<-valid_pathways
        
        for(path in rownames(mat)){
            mat[path,]<-r$misc[['pwpca']][[path]]$xp$scores
        }
        matrices_pagoda2[[names[i]]]<-mat
        
        i=i+1
    }
}

for(name in c("kidney_hallmark","crc_kegg","crc_reactome")){
    write.table(matrices_pagoda2[[name]],file=paste0(outdir_files,name,"_pagoda2.tsv"),col.names = T,row.names = T,quote=F)
}


## AUCell

matrices_aucell<-list()

# kidney
set.seed(123)
cells_rankings <- AUCell_buildRankings(logcpm_kidney)

## hallmark
cells_AUC <- AUCell_calcAUC(hallmark_modules, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

cells_AUC<-getAUC(cells_AUC)
matrices_aucell[["kidney_hallmark_full"]]<-cells_AUC

# CRC
set.seed(123)
cells_rankings <- AUCell_buildRankings(logcpm_crc)

## kegg
cells_AUC <- AUCell_calcAUC(kegg_modules, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

cells_AUC<-getAUC(cells_AUC)
matrices_aucell[["crc_kegg_full"]]<-cells_AUC

## Reactome
cells_AUC <- AUCell_calcAUC(reactome_modules, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

cells_AUC<-getAUC(cells_AUC)
matrices_aucell[["crc_reactome_full"]]<-cells_AUC

# save

for(name in names(matrices_aucell)){
    write.table(matrices_aucell[[name]],file=paste0(outdir_files,name,"_aucell.tsv"),col.names = T,row.names = T,quote=F)
}


## Dilution experiment ##

# Allograft rejection
# baseline

ref_obj<-MAYA_pathway_analysis(expr_mat = logcpm_kidney,
                               modules_list = hallmark_modules["HALLMARK_ALLOGRAFT_REJECTION"],
                               min_cells_pct = 0.05,
                               is_logcpm = F,
                               max_contrib=0.5,
                               compute_umap=F,
                               scale_before_pca = T,
                               all_PCs_in_range = F)

ref_mat<-build_activity_mat(PCA_object = ref_obj$PCA_obj,scaled = T)
ref_spe<-specificity_table(ref_mat,meta_kidney,"Cell_type")
saveRDS(ref_spe,file=file.path(outdir_files,"ref_spe_allograft.rds"))


# Build gene list
{
    # initial pathway genes
    allograft_genes<-hallmark_modules[["HALLMARK_ALLOGRAFT_REJECTION"]]
    
    # genes we draw random genes from
    available_genes<-setdiff(rownames(logcpm_kidney[rowSums(logcpm_kidney!=0)>=10, ]),allograft_genes)
    
    i=1
    seeds<-1:100
    noise10<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise10[[paste0(k,"_AR")]]<-c(allograft_genes,
                                      sample(x=available_genes,size = 10,replace = F))
        i=i+1
        
    }
    names(noise10)<-paste0("noise10_",names(noise10))
    
    i=1
    seeds<-1:100
    noise50<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise50[[paste0(k,"_AR")]]<-c(allograft_genes,
                                      sample(x=available_genes,size = 50,replace = F))
        i=i+1
        
    }
    names(noise50)<-paste0("noise50_",names(noise50))
    
    i=1
    seeds<-1:100
    noise100<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise100[[paste0(k,"_AR")]]<-c(allograft_genes,
                                       sample(x=available_genes,size = 100,replace = F))
        i=i+1
        
    }
    names(noise100)<-paste0("noise100_",names(noise100))
    
    i=1
    seeds<-1:100
    noise200<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise200[[paste0(k,"_AR")]]<-c(allograft_genes,
                                       sample(x=available_genes,size = 200,replace = F))
        i=i+1
        
    }
    names(noise200)<-paste0("noise200_",names(noise200))
    
}

PCA_obj<-MAYA_pathway_analysis(expr_mat = logcpm_kidney,
                                     modules_list = c(noise10,noise50,noise100,noise200),
                                     min_cells_pct = 0.05,
                                     is_logcpm = F,
                                     max_contrib=0.5,
                                     compute_umap=F,
                               scale_before_pca = T,
                               all_PCs_in_range = F)

saveRDS(object = PCA_obj,file=file.path(outdir_files,"random_allograft.rds"))

# Save dataframes for plots

# boxplot
{
    # Build specificity table for each mode
    kidney_noise_maya<-build_activity_mat(PCA_object = PCA_obj$PCA_obj,scaled = T)
    spe_maya<-specificity_table(kidney_noise_maya,meta_kidney,"Cell_type")
    spe_maya<-as.matrix(spe_maya[order(rownames(spe_maya)),])
    
    spe_maya<-spe_maya[order(rownames(spe_maya)),]
    
    PC1_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"_mode")[[1]][2])=="1")]
    PC2_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="2")]
    
    order(c("noise10","noise50","noise100","noise200"))
    PC1_noise_labels<-sapply(PC1_modes,function(x) strsplit(x,"_")[[1]][1])
    PC1_avg<-spe_maya[PC1_modes,"S_MNP1"]
    
    PC2_noise_labels<-sapply(PC2_modes,function(x) strsplit(x,"_")[[1]][1])
    PC2_avg<-spe_maya[PC2_modes,"S_TCD8"]
    
    df<-data.frame(cell_type=c(rep("MNP1",length(PC1_noise_labels)),rep("TCD8",length(PC2_noise_labels))),
                   noise_pct=c(PC1_noise_labels,PC2_noise_labels),
                   maya=c(PC1_avg,PC2_avg))
    df$noise_pct <- factor(df$noise_pct , levels=c("noise10", "noise50", "noise100", "noise200"))
    saveRDS(df,file=file.path(outdir_files,"df_allograft.rds"))
}

# barplot
{
    spe_maya<-specificity_table(kidney_noise_maya,meta_kidney,"Cell_type")
    spe_maya<-as.matrix(spe_maya[order(rownames(spe_maya)),])
    
    PC1_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"_mode")[[1]][2])=="1")]
    PC2_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="2")]
    PC3_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="3")]
    
    noise<-c(sapply(PC1_modes,function(x) strsplit(x,"_")[[1]][1]),
             sapply(PC2_modes,function(x) strsplit(x,"_")[[1]][1]),
             sapply(PC3_modes,function(x) strsplit(x,"_")[[1]][1]))
    
    colnames(spe_maya)<-sapply(colnames(spe_maya),function(x) strsplit(x,"S_")[[1]][2])
    spe_maya<-spe_maya[names(noise),]
    cell_type<-c(colnames(spe_maya)[apply(spe_maya[PC1_modes,],1,which.max)],
                 colnames(spe_maya)[apply(spe_maya[PC2_modes,],1,which.max)],
                 colnames(spe_maya)[apply(spe_maya[PC3_modes,],1,which.max)])
    
    df<-data.frame(noise=as.vector(unlist(noise)),cell_type=as.vector(cell_type),mode=paste0("mode",as.vector(unlist(sapply(names(noise),function(x) strsplit(x,"_mode")[[1]][2])))))
    final <-df %>% dplyr::count(noise, cell_type,mode)
    final$noise <- ordered(final$noise, levels = c("noise10","noise50","noise100","noise200"))
    saveRDS(final,file=file.path(outdir_files,"final_noise_allograft.rds"))
}



# TNFA signaling

# Build gene list
{
    # initial pathway genes
    allograft_genes<-hallmark_modules[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]
    
    # genes we draw random genes from
    available_genes<-setdiff(rownames(logcpm_kidney[rowSums(logcpm_kidney!=0)>=10, ]),allograft_genes)
    
    i=1
    seeds<-1:100
    noise10<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise10[[paste0(k,"_AR")]]<-c(allograft_genes,
                                      sample(x=available_genes,size = 10,replace = F))
        i=i+1
        
    }
    names(noise10)<-paste0("noise10_",names(noise10))
    
    i=1
    seeds<-1:100
    noise50<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise50[[paste0(k,"_AR")]]<-c(allograft_genes,
                                      sample(x=available_genes,size = 50,replace = F))
        i=i+1
        
    }
    names(noise50)<-paste0("noise50_",names(noise50))
    
    i=1
    seeds<-1:100
    noise100<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise100[[paste0(k,"_AR")]]<-c(allograft_genes,
                                       sample(x=available_genes,size = 100,replace = F))
        i=i+1
        
    }
    names(noise100)<-paste0("noise100_",names(noise100))
    
    i=1
    seeds<-1:100
    noise200<-list()
    for(k in 1:100){
        set.seed(seeds[i])
        noise200[[paste0(k,"_AR")]]<-c(allograft_genes,
                                       sample(x=available_genes,size = 200,replace = F))
        i=i+1
        
    }
    names(noise200)<-paste0("noise200_",names(noise200))
    
}

PCA_obj<-MAYA_pathway_analysis(expr_mat = logcpm_kidney,
                                     modules_list = c(noise10,noise50,noise100,noise200),
                                     min_cells_pct = 0.05,
                                     is_logcpm = F,
                                     max_contrib=0.5,
                                     compute_umap=F,
                               scale_before_pca = T,
                               all_PCs_in_range = F)

saveRDS(object = PCA_obj,file=file.path(outdir_files,"random_TNFA.rds"))

# Save dataframes for plots

# barplot
{
    kidney_noise_maya<-build_activity_mat(PCA_object = PCA_obj$PCA_obj,scaled = T)
    spe_maya<-specificity_table(kidney_noise_maya,meta_kidney,"Cell_type")
    spe_maya<-as.matrix(spe_maya[order(rownames(spe_maya)),])
    
    PC1_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"_mode")[[1]][2])=="1")]
    PC2_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="2")]
    PC3_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="3")]
    PC4_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="4")]
    PC5_modes<-rownames(spe_maya)[which(sapply(rownames(spe_maya),function(x) strsplit(x,"mode")[[1]][2])=="5")]
    
    noise<-c(sapply(PC1_modes,function(x) strsplit(x,"_")[[1]][1]),
             sapply(PC2_modes,function(x) strsplit(x,"_")[[1]][1]),
             sapply(PC3_modes,function(x) strsplit(x,"_")[[1]][1]),
             sapply(PC4_modes,function(x) strsplit(x,"_")[[1]][1]),
             sapply(PC5_modes,function(x) strsplit(x,"_")[[1]][1]))
    
    colnames(spe_maya)<-sapply(colnames(spe_maya),function(x) strsplit(x,"S_")[[1]][2])
    spe_maya<-spe_maya[names(noise),]
    cell_type<-c(colnames(spe_maya)[apply(spe_maya[PC1_modes,],1,which.max)],
                 colnames(spe_maya)[apply(spe_maya[PC2_modes,],1,which.max)],
                 colnames(spe_maya)[apply(spe_maya[PC3_modes,],1,which.max)],
                 colnames(spe_maya)[apply(spe_maya[PC4_modes,],1,which.max)],
                 colnames(spe_maya)[apply(spe_maya[PC5_modes,],1,which.max)])
    
    df<-data.frame(noise=unlist(noise),cell_type,mode=paste0("mode",sapply(names(noise),function(x) strsplit(x,"_mode")[[1]][2])))
    final<-df %>% dplyr::count(noise, cell_type,mode)
    final$noise <- ordered(final$noise, levels = c("noise10","noise50","noise100","noise200"))
    saveRDS(final,file=file.path(outdir_files,"final_noise_TNFA.rds"))
}


#############


### SS2 vs 10X ###


data("pbmcsca")
seurat <- pbmcsca
rm(pbmcsca)
seurat$percent.mito<-as.numeric(seurat$percent.mito)

# QC and subset
seurat<-subset(seurat, subset= nFeature_RNA>200 & percent.mito<20)
seurat_small<-subset(seurat, subset= Method %in% c("Smart-seq2","10x Chromium (v3)"))
saveRDS(seurat_small,file=file.path(outdir_files,"pbmcsca_seurat.rds"))

seurat_SS2<-subset(seurat, subset= Method =="Smart-seq2") # 526 cells
seurat_10X<-subset(seurat, subset= Method =="10x Chromium (v3)") # 2993 cells

# MAYA

SS2_hallmark<-MAYA_pathway_analysis(expr_mat = GetAssayData(seurat_SS2),
                                    modules_list = "hallmark",
                                    min_cells_pct = 0.05,
                                    is_logcpm = F,
                                    scale_before_pca = T,
                                    all_PCs_in_range = F)
saveRDS(SS2_hallmark,file=file.path(outdir_files,"pbmcsca_SS2_hallmark.rds"))

tenX_hallmark<-MAYA_pathway_analysis(expr_mat = GetAssayData(seurat_10X),
                                     modules_list = "hallmark",
                                     min_cells_pct = 0.05,
                                     is_logcpm = F,
                                     scale_before_pca = T,
                                     all_PCs_in_range = F)
saveRDS(tenX_hallmark,file=file.path(outdir_files,"pbmcsca_tenX_hallmark.rds"))



######################


rm(list=ls())
gc()























