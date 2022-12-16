### Plots Figure 2 ###

indir_files="./output_dir/"
outdir_plots="./output_plots/Figure2/"
dir.create(outdir_plots)

source("./scripts/utils.R")

#################################


## load result file 
kidney_hallmark<-readRDS(file.path(indir_files,"PCA_obj_hallmark_kidney.rds"))
meta_kidney<-readRDS(file.path(indir_files,"meta_kidney.rds"))
logcpm_kidney<-readRDS(file.path(indir_files,"logcpm_kidney.rds"))

pagoda2<-as.matrix(read.table(file.path(indir_files,"kidney_hallmark_pagoda2.tsv"),header = T,check.names = F))
aucell<-as.matrix(read.table(file.path(indir_files,"kidney_hallmark_full_aucell.tsv"),header = T,check.names = F))

activity_mat<-build_activity_mat(kidney_hallmark$PCA_obj,scaled = T)

### Plot ###
    
# 2.a Activity matrix kidney hallmark scaled

p<-plot_heatmap_activity_mat(activity_mat,meta = meta_kidney,annot_name = "Cell_type")
save_plot(p,file.path(outdir_plots,"2.A.kidney_heatmap_hallmark_activity_matrix_scaled.png"),width = 15,height = 8)

# 2.b PCA plot mode 1 vs mode 2

module="HALLMARK_ALLOGRAFT_REJECTION"
PCA=kidney_hallmark$PCA_obj[[module]]$activity_scores_raw[1:2,]

df<-data.frame(PC1=PCA[1,],PC2=PCA[2,],Cell_type=meta_kidney[,"Cell_type"])
p <- ggplot(df, aes(x=PC1, y=PC2,color=Cell_type)) + geom_point(shape=21,colour="white",aes(fill=Cell_type),size=2,stroke=0.5) + theme_classic() + scale_fill_manual(values=generate_palette(5),breaks=unique(meta_kidney[,"Cell_type"]))+ theme(legend.position = "none",plot.title=element_blank())
p <- ggExtra::ggMarginal(p, type = "histogram",xparams = list(binwidth = 1))
save_plot(p,file.path(outdir_plots,"2.B.kidney_PCA_hallmark_allo_rej_PC12_no_legend.png"),width = 7,height = 7)


# 2.c top contributing genes allograft rejection

module="HALLMARK_ALLOGRAFT_REJECTION"
p<-plot_top_contrib(expr_mat=logcpm_kidney, PCA_object=kidney_hallmark$PCA_obj, module=module, n = 10, meta = meta_kidney,annot_name = "Cell_type")
save_plot(p,file.path(outdir_plots,"2.C.kidney_heatmap_hallmark_allo_rej_top_genes.png"),width = 15,height = 11)


# 2.d UMAP MAYA colored by annotation and by allograft rejection activity

p<-plot_umap_annot(kidney_hallmark$umap,type="Cell type",labels=meta_kidney[,"Cell_type"]) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"2.D.kidney_UMAP_author_annot_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(kidney_hallmark$umap,expr_mat = activity_mat,gene = "HALLMARK_ALLOGRAFT_REJECTION_mode1") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"2.D.kidney_UMAP_allo_rej_mode1_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(kidney_hallmark$umap,expr_mat = activity_mat,gene = "HALLMARK_ALLOGRAFT_REJECTION_mode2") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"2.D.kidney_UMAP_allo_rej_mode2_no_legend.png"),width = 7,height = 7)


# 2.e Heatmap comparing MAYA modes to AUCell and Pagoda2 for allograft rejection

mat=matrix(c(kidney_hallmark$activity_matrix["HALLMARK_ALLOGRAFT_REJECTION_mode1",],
             kidney_hallmark$activity_matrix["HALLMARK_ALLOGRAFT_REJECTION_mode2",],
             aucell["HALLMARK_ALLOGRAFT_REJECTION",],
             pagoda2["HALLMARK_ALLOGRAFT_REJECTION",]),nrow=4,ncol=ncol(logcpm_kidney),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2","AUCell","Pagoda2")
colnames(mat)<-colnames(logcpm_kidney)

# Data frame with column annotations.
mat_col <- data.frame(Cell_type = meta_kidney$Cell_type)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(5))
names(mat_colors$Cell_type) <- unique(meta_kidney$Cell_type)

p<-pheatmap(
    mat               = mat,
    color             = viridis(10),
    border_color      = NA,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 7,
    main              = "Kidney - HALLMARK_ALLOGRAFT_REJECTION",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8,
    gaps_col = c(305,305+259,155+305+259,415+305+155+259),gaps_row = c(2,3)
)

save_plot(p,file.path(outdir_plots,"2.E.kidney_heatmap_activity_hallmark_allo_rej.png"),width = 15,height = 8)


# 2.f Noise dilution barplot allograft

final_noise_allograft<-readRDS(file=file.path(indir_files,"final_noise_allograft.rds"))

p<-ggplot(final_noise_allograft, aes(x=noise, y=n, fill=cell_type)) +
       geom_bar(stat = "identity")+
       xlab("noise")+ylab("Number of detected mode")+
       ggtitle("Detected modes - ALLOGRAFT REJECTION")+theme_classic()+
       scale_fill_manual(values = c(generate_palette(5)),breaks=unique(meta_kidney$Cell_type))+
       facet_wrap(~mode)+ theme(plot.title = element_text(hjust = 0.5))+
       theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))
save_plot(p,file.path(outdir_plots,"2.F.noise_allo_rej_summary.png"),width = 15,height = 11)

# 2.g Noise dilution boxplot

df_allograft<-readRDS(file=file.path(indir_files,"df_allograft.rds"))

p<-ggboxplot(
    as_tibble(df_allograft),
    x = "noise_pct",
    y= "maya",
    width = 0.5
) + geom_point(aes(fill = cell_type, x = noise_pct, y = maya), shape = 21, size = 4,position = position_jitterdodge(jitter.width = 0.3))+
    facet_wrap(~cell_type)+theme_classic()+
    theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8)
    ) +scale_fill_manual(values=c(generate_palette(4)[4],generate_palette(4)[3]),breaks=c("MNP1","TCD8"))+
    ggtitle("Specicity score by cell type with noise")+
    geom_hline(data=data.frame(yint=0.566,cell_type="MNP1"),aes(yintercept=yint), linetype="dashed", color = "black")+
    geom_hline(data=data.frame(yint=0.879,cell_type="TCD8"),aes(yintercept=yint), linetype="dashed", color = "black")

save_plot(p,file.path(outdir_plots,"2.G.noise_allo_rejection_spe_score_no_boxplot.png"),width = 8,height = 8)







# Supp

outdir_plots="./output_plots/Figure2_sup/"
dir.create(outdir_plots)

# sup2.a PCA plot mode 1,2,3,4 TNFA pathway 

module="HALLMARK_TNFA_SIGNALING_VIA_NFKB"
PCA=kidney_hallmark$PCA_obj[[module]]$activity_scores_raw[1:4,]

df<-data.frame(PC1=PCA[1,],PC2=PCA[2,],Cell_type=meta_kidney[,"Cell_type"])
p <- ggplot(df, aes(x=PC1, y=PC2,color=Cell_type)) + geom_point(shape=21,colour="white",aes(fill=Cell_type),size=2,stroke=0.5) + theme_classic() + scale_fill_manual(values=generate_palette(5),breaks=unique(meta_kidney[,"Cell_type"]))+ theme(legend.position = "none",plot.title=element_blank())
p <- ggExtra::ggMarginal(p, type = "histogram",xparams = list(binwidth = 1))
save_plot(p,file.path(outdir_plots,"sup2.A.kidney_PCA_hallmark_TNFA_PC12_no_legend.png"),width = 7,height = 7)


df<-data.frame(PC3=PCA[3,],PC4=PCA[4,],Cell_type=meta_kidney[,"Cell_type"])
p <- ggplot(df, aes(x=PC3, y=PC4,color=Cell_type)) + geom_point(shape=21,colour="white",aes(fill=Cell_type),size=2,stroke=0.5) + theme_classic() + scale_fill_manual(values=generate_palette(5),breaks=unique(meta_kidney[,"Cell_type"]))+ theme(legend.position = "none",plot.title=element_blank())
p <- ggExtra::ggMarginal(p, type = "histogram",xparams = list(binwidth = 1))
save_plot(p,file.path(outdir_plots,"sup2.A.kidney_PCA_hallmark_TNFA_PC34_no_legend.png"),width = 7,height = 7)


# sup2.b heatmap top contributing genes TNFA

module="HALLMARK_TNFA_SIGNALING_VIA_NFKB"
p<-plot_top_contrib(expr_mat=logcpm_kidney, PCA_object=kidney_hallmark$PCA_obj, module=module, n = 10, meta = meta_kidney,annot_name = "Cell_type")
save_plot(p,file.path(outdir_plots,"sup2.B.kidney_heatmap_hallmark_allo_rej_top_genes.png"),width = 15,height = 11)


# sup2.c UMAP colored by TNFA activity four modes

p<-plot_umap_gene(kidney_hallmark$umap,expr_mat = activity_mat,gene = "HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode1") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup2.C.kidney_UMAP_TNFA_mode1_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(kidney_hallmark$umap,expr_mat = activity_mat,gene = "HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode2") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup2.C.kidney_UMAP_TNFA_mode2_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(kidney_hallmark$umap,expr_mat = activity_mat,gene = "HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode3") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup2.C.kidney_UMAP_TNFA_mode3_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(kidney_hallmark$umap,expr_mat = activity_mat,gene = "HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode4") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup2.C.kidney_UMAP_TNFA_mode4_no_legend.png"),width = 7,height = 7)


# sup2.d heatmap comparing MAYA to AUCell and Pagoda2 for TNFA

mat=matrix(c(kidney_hallmark$activity_matrix["HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode1",],
             kidney_hallmark$activity_matrix["HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode2",],
             kidney_hallmark$activity_matrix["HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode3",],
             kidney_hallmark$activity_matrix["HALLMARK_TNFA_SIGNALING_VIA_NFKB_mode4",],
             aucell["HALLMARK_TNFA_SIGNALING_VIA_NFKB",],
             pagoda2["HALLMARK_TNFA_SIGNALING_VIA_NFKB",]),nrow=6,ncol=ncol(logcpm_kidney),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2","MAYA_mode_3","MAYA_mode_4","AUCell","Pagoda2")
colnames(mat)<-colnames(logcpm_kidney)

# Data frame with column annotations.
mat_col <- data.frame(Cell_type = meta_kidney$Cell_type)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(5))
names(mat_colors$Cell_type) <- unique(meta_kidney$Cell_type)

p<-pheatmap(
    mat               = mat,
    color             = viridis(10),
    border_color      = NA,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 7,
    main              = "Kidney - HALLMARK_TNFA_SIGNALING_VIA_NFKB ",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8,
    gaps_col = c(305,305+259,155+305+259,415+305+155+259),gaps_row = c(4,5)
)
save_plot(p,file.path(outdir_plots,"sup2.D.kidney_heatmap_activity_TNFA_rej.png"),width = 15,height = 8)


# sup2.e Noise dilution barplot

final_noise_TNFA<-readRDS(file=file.path(indir_files,"final_noise_TNFA.rds"))

p<-ggplot(final_noise_TNFA, aes(x=noise, y=n, fill=cell_type)) +
    geom_bar(stat = "identity")+
    xlab("noise")+ylab("Number of detected mode")+
    ggtitle("Detected modes - TNFA SIGNALING VIA NFKB")+theme_classic()+
    scale_fill_manual(values = c(generate_palette(5)),breaks=unique(meta_kidney$Cell_type))+
    facet_wrap(~mode,ncol=5)+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5))

save_plot(p,file.path(outdir_plots,"sup2.E.noise_TNFA_summary.png"),width = 15,height = 11)




rm(list=ls())
gc()









