### Plots Figure 5 ###

indir_files="./output_dir/"
outdir_plots="./output_plots/Figure5/"
dir.create(outdir_plots)

source("./scripts/utils.R")

### Load data ###

# hallmark
umap_ovary_gene<-readRDS(file.path(indir_files,"umap_ovary_gene_based.rds"))
logcpm_ovary<-readRDS(file.path(indir_files,"logcpm_ovary.rds"))
ovary_hallmark<-readRDS(file.path(indir_files,"PCA_obj_hallmark_ovary.rds"))
meta_ovary<-readRDS(file.path(indir_files,"meta_ovary.rds"))
final_df_patient_ovary<-read.table(file=paste0(indir_files,"ovary_shannon_results.tsv"),sep="\t",header=T)

act_mat=build_activity_mat(ovary_hallmark$PCA_obj,scaled=T)
avg<-average_by_cluster(act_mat,meta_ovary[colnames(act_mat),"Cell_type"])

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", ""))


order_cells<-data.frame(cells=rownames(meta_ovary),subtype=meta_ovary$Cell_type)
order_cells <- order_cells %>%
    mutate(category =  factor(subtype, levels = unique(meta_ovary$Cell_type))) %>%
    arrange(category)


meta_ovary_tmp<-meta_ovary[order_cells$cells,]
logcpm_ovary_tmp<-logcpm_ovary[,order_cells$cells]

PCA_obj<-ovary_hallmark$PCA_obj[grep("_DN$",x = names(ovary_hallmark$PCA_obj),invert = T,value = T)]

# KEGG
ovary_kegg<-readRDS(file=file.path(indir_files,"ovary_kegg.rds"))

# GSEA


# NMF
nmf_intersect_melt_tumor<-readRDS(file.path(indir_files,"NMF_correlation_mat.rds"))
matrix_pct<-readRDS(file.path(outdir_files,"NMF_percentage_with_modes.rds"))
NMF_enrichment<-readRDS(file.path(indir_files,"NMF_enrichment.rds"))

# other datasets

meta_breast<-readRDS(file.path(indir_files,"meta_breast.rds"))
logcpm_breast<-readRDS(file.path(indir_files,"logcpm_breast.rds"))
breast_hallmark<-readRDS(file.path(indir_files,"breast_hallmark.rds"))

meta_lung<-readRDS(file.path(indir_files,"meta_lung.rds"))
logcpm_lung<-readRDS(file.path(indir_files,"logcpm_lung.rds"))
lung_hallmark<-readRDS(file.path(indir_files,"lung_hallmark.rds"))



###



# 5.a UMAP MAYA based author annot and patient
p<-plot_umap_annot_random(ovary_hallmark$umap,labels = meta_ovary$Cell_type) + theme(legend.position = "none",plot.title =element_blank())
save_plot(p,file.path(outdir_plots,"5.A.Ovary_UMAP_author_annotation_HALLMARK.png"),width = 7,height = 7)

p<-plot_umap_annot_random(ovary_hallmark$umap,labels = meta_ovary$patient_id,type = "Patient",colors = generate_palette(30)[14:24]) + theme(legend.position = "none",plot.title =element_blank())
save_plot(p,file.path(outdir_plots,"5.A.Ovary_UMAP_patient_id_HALLMARK.png"),width = 7,height = 7)

# 5.b. SDI with tumor clusters highlighted
{
    
    citeCol <- c(wes_palette("IsleofDogs1")[6],wes_palette("GrandBudapest1")[2])
    
    names(citeCol) <- c("FALSE","TRUE")
    
    p<-ggboxplot(
        as_tibble(final_df_patient_ovary),
        x = "method",
        y= "shannon",
        width = 0.5
    ) + geom_point(aes(fill = is_tumor, x = method, y = shannon), shape = 21, size = 5,position = position_jitterdodge())+
        theme(
            aspect.ratio = 1.5,
            legend.title = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 8, face = "bold"),
            axis.text.y = element_text(size = 8, face = "bold"),
            plot.title = element_text(hjust = 0.5)
        ) + grids(linetype = "dashed")+ggtitle("Shannon index by cluster for different methods")+
        scale_fill_manual(values = citeCol)+
        scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.1))+ theme(plot.title =element_blank())+
        stat_compare_means(label = "p.signif", method = "wilcox.test",
                           ref.group = "MAYA",symnum.args=symnum.args,size=3)
    
}
save_plot(p,file.path(outdir_plots,"5.B.Ovary_shannon_index.png"),width = 15,height = 11)



# 5.c Barplot top activated pathway for epithelial, macrophages, CAFs and T cells


summary_celltype_hallmark<-generate_summary_specific_modules(PCA_obj,meta = meta_ovary,annot_name = "Cell_type",top_genes = 0)
generate_summary_specific_modules(PCA_obj,meta = meta_ovary,annot_name = "Cell_type",top_genes = 5,file = file.path(outdir_plots,"summary_top_spe_modes.xlsx"))

tmp<-head(summary_celltype_hallmark[["EOC"]],n=5)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",5))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(1)),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1)) + coord_flip()
save_plot(p,file.path(outdir_plots,"5.C.ovary_barplot_epithelial.png"),width = 7,height = 3)

tmp<-head(summary_celltype_hallmark[["Macrophages"]],n=5)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",5))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(8)[8]),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1))+ coord_flip()
save_plot(p,file.path(outdir_plots,"5.C.ovary_barplot_macrophages.png"),width = 7,height = 3)

tmp<-head(summary_celltype_hallmark[["CAF"]],n=5)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",5))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(2)[2]),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1))+ coord_flip()
save_plot(p,file.path(outdir_plots,"5.C.ovary_barplot_CAF.png"),width = 7,height = 3)


tmp<-head(summary_celltype_hallmark[["T-cells"]],n=5)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",5))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(5)[5]),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1))+ coord_flip()
save_plot(p,file.path(outdir_plots,"5.C.ovary_barplot_Tcells.png"),width = 7,height = 3)


# 5.d Heatmap MAYA modes EMT in all cells

mat=matrix(c(ovary_hallmark$activity_matrix["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode1",],
             ovary_hallmark$activity_matrix["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode2",],
             ovary_hallmark$activity_matrix["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode3",],
             ovary_hallmark$activity_matrix["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode4",],
             ovary_hallmark$activity_matrix["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode5",]),
           nrow=5,ncol=ncol(logcpm_ovary),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2","MAYA_mode_3","MAYA_mode_4","MAYA_mode_5")
colnames(mat)<-colnames(ovary_hallmark$activity_matrix)

mat<-mat[,order_cells$cells]

mat_col <- data.frame(Cell_type = meta_ovary[order_cells$cells,]$Cell_type)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(13))
names(mat_colors$Cell_type) <- unique(meta_ovary$Cell_type)

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
    main              = "Ovary - HALLMARK_EMT",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8
)
save_plot(p,file.path(outdir_plots,"5.D.ovary_heatmap_activity_hallmark_EMT.png"),width = 13,height = 4)

module="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
p<-plot_top_contrib(expr_mat=logcpm_ovary_tmp, PCA_object=ovary_hallmark$PCA_obj, module=module, n = 10, meta = meta_ovary_tmp,annot_name = "Cell_type",cluster_cols = F)
save_plot(p,file.path(outdir_plots,"5.D.ovary_heatmap_hallmark_EMT_top_genes.png"),width = 15,height = 12)



# 5.e UMAP MAYA based with EMT activity for 3 modes and corresponding violin plot with activity by cell type
mode="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode1"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"5.E.UMAP_EMT_mode1.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"5.E.Vln_EMT_mode1_stat.png"),width = 14,height = 8)

mode="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode2"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"5.E.UMAP_EMT_mode2.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"5.E.Vln_EMT_mode2_stat.png"),width = 14,height = 8)

mode="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_mode3"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"5.E.UMAP_EMT_mode3.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"5.E.Vln_EMT_mode3_stat.png"),width = 14,height = 8)


# 5.f heatmap MAYA modes estrogen and top contributors

mat=matrix(c(ovary_hallmark$activity_matrix["HALLMARK_ESTROGEN_RESPONSE_EARLY_mode1",],
             ovary_hallmark$activity_matrix["HALLMARK_ESTROGEN_RESPONSE_EARLY_mode2",]),
           nrow=2,ncol=ncol(logcpm_ovary),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2")
colnames(mat)<-colnames(ovary_hallmark$activity_matrix)

mat<-mat[,order_cells$cells]

# Data frame with column annotations.
mat_col <- data.frame(Cell_type = meta_ovary[order_cells$cells,]$Cell_type)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(13))
names(mat_colors$Cell_type) <- unique(meta_ovary$Cell_type)

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
    main              = "Ovary - HALLMARK_ESTROGEN",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8
)
save_plot(p,file.path(outdir_plots,"5.F.ovary_heatmap_activity_hallmark_estrogen.png"),width = 13,height = 2)


module="HALLMARK_ESTROGEN_RESPONSE_EARLY"
p<-plot_top_contrib(expr_mat=logcpm_ovary_tmp, PCA_object=ovary_hallmark$PCA_obj, module=module, n = 10, meta = meta_ovary_tmp,annot_name = "Cell_type",cluster_cols = F)
save_plot(p,file.path(outdir_plots,"5.F.ovary_heatmap_hallmark_estrogen_top_genes.png"),width = 15,height = 11)


# 5.g heatmap MAYA modes coagulation and top contributors

mat=matrix(c(ovary_hallmark$activity_matrix["HALLMARK_COAGULATION_mode1",],
             ovary_hallmark$activity_matrix["HALLMARK_COAGULATION_mode2",],
             ovary_hallmark$activity_matrix["HALLMARK_COAGULATION_mode3",]),
           nrow=3,ncol=ncol(logcpm_ovary),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2","MAYA_mode_3")
colnames(mat)<-colnames(ovary_hallmark$activity_matrix)

mat<-mat[,order_cells$cells]

# Data frame with column annotations.
mat_col <- data.frame(Cell_type = meta_ovary[order_cells$cells,]$Cell_type)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(13))
names(mat_colors$Cell_type) <- unique(meta_ovary$Cell_type)

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
    main              = "Ovary - HALLMARK_COAG",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8
)
save_plot(p,file.path(outdir_plots,"5.G.ovary_heatmap_activity_hallmark_coagulation.png"),width = 13,height = 3)


module="HALLMARK_COAGULATION"
p<-plot_top_contrib(expr_mat=logcpm_ovary_tmp, PCA_object=ovary_hallmark$PCA_obj, module=module, n = 10, meta = meta_ovary_tmp,annot_name = "Cell_type",cluster_cols = F)
save_plot(p,file.path(outdir_plots,"5.G.ovary_heatmap_hallmark_coagulation_top_genes.png"),width = 15,height = 11)





# 5.f heatmap MAYA modes Kegg WNT and top contributors

mat=matrix(c(ovary_kegg$activity_matrix["KEGG_WNT_SIGNALING_PATHWAY_mode1",],
             ovary_kegg$activity_matrix["KEGG_WNT_SIGNALING_PATHWAY_mode2",]),
           nrow=2,ncol=ncol(logcpm_ovary),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2")
colnames(mat)<-colnames(ovary_kegg$activity_matrix)

mat<-mat[,order_cells$cells]

mat_col <- data.frame(Cell_type = meta_ovary[order_cells$cells,]$Cell_type)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(13))
names(mat_colors$Cell_type) <- unique(meta_ovary$Cell_type)

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
    main              = "Ovary - KEGG WNT",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8
)
save_plot(p,file.path(outdir_plots,"5.H.ovary_kegg_wnt_activity.png"),width = 13,height = 2)


module="KEGG_WNT_SIGNALING_PATHWAY"
p<-plot_top_contrib(expr_mat=logcpm_ovary_tmp, PCA_object=ovary_kegg$PCA_obj, module=module, n = 10, meta = meta_ovary_tmp,annot_name = "Cell_type",cluster_cols = F)
save_plot(p,file.path(outdir_plots,"5.H.ovary_kegg_wnt_top_contrib.png"),width = 15,height = 11)









# Supp


outdir_plots="./output_plots/Figure5_sup/"
dir.create(outdir_plots)


# sup5.a Activity matrix hallmark ovary

{
    meta=meta_ovary
    annot_name=c("Cell_type","patient_id")
    clustering_distance="euclidean"
    clustering_method="ward.D2"
    fontsize_row=5
    fontsize_col=5
    activity_mat=build_activity_mat(ovary_hallmark$PCA_obj)
    
    spe_table<-specificity_table(build_activity_mat(ovary_hallmark$PCA_obj),meta = meta_ovary,annot_name = "Cell_type")
    maxs<-rowMax(as.matrix(spe_table))
    activity_mat<-activity_mat[maxs>(1/13)*1.5,]
    
    tmp<-sapply(rownames(activity_mat),function(x) gsub("HALLMARK_","",x))
    rownames(activity_mat)=tmp
    
    mat_col=meta[,annot_name]
    
    mat_colors <- list()
    mat_colors[[annot_name[1]]]<-generate_palette(length(unique(meta[,annot_name[1]])))
    names(mat_colors[[annot_name[1]]]) <- unique(meta[,annot_name[1]])
    mat_colors[[annot_name[2]]]<-generate_palette(length(unique(meta[,annot_name[1]]))+length(unique(meta[,annot_name[2]])))[(length(unique(meta[,annot_name[1]]))+1):(length(unique(meta[,annot_name[1]]))+length(unique(meta[,annot_name[2]])))]
    names(mat_colors[[annot_name[2]]]) <- unique(meta[,annot_name[2]])
    
    p<-pheatmap(
        activity_mat,
        annotation_col = mat_col,
        cluster_cols = T,
        cluster_rows=T,
        fontsize = 5,
        clustering_distance_rows = clustering_distance,
        clustering_distance_cols = clustering_distance,
        clustering_method=clustering_method,
        show_colnames = F,
        annotation_colors = mat_colors,
        color             = viridis(10,direction=1),
        main= "Activity matrix",
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col)
}
save_plot(p,file.path(outdir_plots,"sup5.A.ovary_heatmap_hallmark_activity_matrix_scaled_filtered.png"),width = 20,height = 15)


# sup5.b UMAP gene based author annot and patient

p<-plot_umap_annot_random(umap_ovary_gene,labels = meta_ovary$Cell_type) + theme(legend.position = "none",plot.title =element_blank())
save_plot(p,file.path(outdir_plots,"sup5.B.Ovary_UMAP_author_annotation_gene_based.png"),width = 7,height = 7)

p<-plot_umap_annot_random(umap_ovary_gene,labels = meta_ovary$patient_id,type = "Patient",colors = generate_palette(30)[14:24]) + theme(legend.position = "none",plot.title =element_blank())
save_plot(p,file.path(outdir_plots,"sup5.B.Ovary_UMAP_patient_id_gene_based.png"),width = 7,height = 7)


# sup5.c comparison GSEA/MAYA

colors<-c(generate_palette(1),generate_palette(2)[2],generate_palette(8)[8])
names(colors)<-c("EOC","CAF","Macrophages")

cell_type<-"EOC"
gsea_res<-readRDS(file.path(indir_files,paste0("gsea_res_",cell_type,".rds")))
gsea_res<-gsea_res[which(gsea_res$up_down=="Up"),]

p<-ggplot(data=gsea_res, aes(x=reorder(pathway, NES), y=logqval, fill=up_down)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(colors[[cell_type]],"grey"),breaks=c("Up","Down"))+
    theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())+
    coord_flip()+ylim(c(0,16))
save_plot(p,file.path(outdir_plots,paste0("sup5.C.",cell_type,"_GSEA_plot_no_legend.png")),width = 15,height = 2)


cell_type<-"CAF"
gsea_res<-readRDS(file.path(indir_files,paste0("gsea_res_",cell_type,".rds")))
gsea_res<-gsea_res[which(gsea_res$up_down=="Up"),]

p<-ggplot(data=gsea_res, aes(x=reorder(pathway, NES), y=logqval, fill=up_down)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(colors[[cell_type]],"grey"),breaks=c("Up","Down"))+
    theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())+
    coord_flip()+ylim(c(0,16))
save_plot(p,file.path(outdir_plots,paste0("sup5.C.",cell_type,"_GSEA_plot_no_legend.png")),width = 15,height = 4)


cell_type<-"Macrophages"
gsea_res<-readRDS(file.path(indir_files,paste0("gsea_res_",cell_type,".rds")))
gsea_res<-gsea_res[which(gsea_res$up_down=="Up"),]

p<-ggplot(data=gsea_res, aes(x=reorder(pathway, NES), y=logqval, fill=up_down)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(colors[[cell_type]],"grey"),breaks=c("Up","Down"))+
    theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())+
    coord_flip()+ylim(c(0,16))
save_plot(p,file.path(outdir_plots,paste0("sup5.C.",cell_type,"_GSEA_plot_no_legend.png")),width = 15,height = 4)


# Compare with MAYA

summary_celltype_hallmark<-generate_summary_specific_modules(PCA_obj,meta = meta_ovary,annot_name = "Cell_type",top_genes = 0)

n=12
tmp<-head(summary_celltype_hallmark[["EOC"]],n=n)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",n))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(1)),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1)) + coord_flip()
save_plot(p,file.path(outdir_plots,"sup5.C.ovary_barplot_epithelial_no_legend.png"),width = 4,height = 3)

n=8
tmp<-head(summary_celltype_hallmark[["Macrophages"]],n=n)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",n))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(8)[8]),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1)) + coord_flip()
save_plot(p,file.path(outdir_plots,"sup5.C.ovary_barplot_macrophages_no_legend.png"),width = 4,height = 2.5)


n=5
tmp<-head(summary_celltype_hallmark[["CAF"]],n=n)
df<-data.frame(mode=sapply(names(tmp),function(x) gsub("HALLMARK_","",x)),specificity=tmp,cell_type=rep("type",n))
p<-ggplot(data=df, aes(x=reorder(mode, specificity), y=specificity, fill=cell_type)) +
    geom_bar(stat="identity",width=0.8)+ theme_minimal()+
    theme_classic()+scale_fill_manual(values = c(generate_palette(2)[2]),breaks="type")+
    theme(legend.position = "none",axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
    ylim(c(0,1)) + coord_flip()
save_plot(p,file.path(outdir_plots,"sup5.C.ovary_barplot_caf_no_legend.png"),width = 4,height = 2)



# sup5.d EMT modes in ovary, breast and lung dataset and overlap between contributing genes

# breast
meta_breast$Annotation<-"Immune other"
meta_breast$Annotation[which(meta_breast$cell_type %in% c("Malignant"))]<-"Tumor cells"
meta_breast$Annotation[which(meta_breast$cell_type %in% c("Fibroblast"))]<-"CAF"
meta_breast$Annotation[which(meta_breast$cell_type %in% c("Macrophage"))]<-"Macrophages"
meta_breast$Annotation[which(meta_breast$cell_type %in% c("Endothelial"))]<-"Endothelial"

order_cells<-data.frame(cells=rownames(meta_breast),subtype=meta_breast$Annotation)
order_cells <- order_cells %>% 
    mutate(category =  factor(subtype, levels = c("Tumor cells","CAF","Macrophages","Endothelial","Immune others"))) %>%
    arrange(category)

module="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
p<-plot_heatmap_activity_mat(build_activity_mat(breast_hallmark$PCA_obj[module],scaled = T)[,order_cells$cells],meta = meta_breast[order_cells$cells,],annot_name = "Annotation",cluster_cols = F,cluster_rows = F)
save_plot(p,file.path(outdir_plots,"sup5.D.breast_activity_EMT.png"),width = 15,height = 2)

# lung
meta_lung$Annotation<-"Immune other"
meta_lung$Annotation[which(meta_lung$cell_type %in% c("Malignant","Epithelial"))]<-"Tumor cells"
meta_lung$Annotation[which(meta_lung$cell_type %in% c("Fibroblast"))]<-"CAF"
meta_lung$Annotation[which(meta_lung$cell_type %in% c("Macrophage"))]<-"Macrophages"
meta_lung$Annotation[which(meta_lung$cell_type %in% c("Endothelial"))]<-"Endothelial"

order_cells<-data.frame(cells=rownames(meta_lung),subtype=meta_lung$Annotation)
order_cells <- order_cells %>% 
    mutate(category =  factor(subtype, levels = c("Tumor cells","CAF","Macrophages","Endothelial","Immune others"))) %>%
    arrange(category)

module="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
p<-plot_heatmap_activity_mat(build_activity_mat(lung_hallmark$PCA_obj[module],scaled = T)[,order_cells$cells],meta = meta_lung[order_cells$cells,],annot_name = "Annotation",cluster_cols = F,cluster_rows = F)
save_plot(p,file.path(outdir_plots,"sup5.D.lung_activity_EMT.png"),width = 15,height = 2)


# ovary
meta_ovary$Annotation<-"Immune other"
meta_ovary$Annotation[which(meta_ovary$Cell_type %in% c("EOC"))]<-"Tumor cells"
meta_ovary$Annotation[which(meta_ovary$Cell_type %in% c("CAF","Mesothelial"))]<-"CAF"
meta_ovary$Annotation[which(meta_ovary$Cell_type %in% c("Macrophages"))]<-"Macrophages"
meta_ovary$Annotation[which(meta_ovary$Cell_type %in% c("Endothelial"))]<-"Endothelial"

order_cells<-data.frame(cells=rownames(meta_ovary),subtype=meta_ovary$Annotation)
order_cells <- order_cells %>% 
    mutate(category =  factor(subtype, levels = c("Tumor cells","CAF","Macrophages","Endothelial","Immune others"))) %>%
    arrange(category)

module="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
p<-plot_heatmap_activity_mat(build_activity_mat(ovary_hallmark$PCA_obj[module],scaled = T)[,order_cells$cells],meta = meta_ovary[order_cells$cells,],annot_name = "Annotation",cluster_cols = F,cluster_rows = F)
save_plot(p,file.path(outdir_plots,"sup5.D.ovary_activity_EMT.png"),width = 15,height = 4)



# EMT CAF
module="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
n=30
top_ovary<-get_top_contrib(ovary_hallmark$PCA_obj,module = module,n=n)
top_breast<-get_top_contrib(breast_hallmark$PCA_obj,module = module,n=n)
top_lung<-get_top_contrib(lung_hallmark$PCA_obj,module = module,n=n)

genes<-unique(unlist(c(top_ovary$PC1,top_breast$PC1,top_lung$PC1)))
mat<-cbind(genes %in% top_ovary$PC1,genes %in% top_breast$PC1,genes %in% top_lung$PC1)
rownames(mat)<-genes
colnames(mat)<-c("ovary","breast","lung")

row_sum<-rowSums(mat)
mat<-mat[order(row_sum,decreasing = T),]
mat[mat ==T]<-1
mat[mat ==F]<-0
p<-pheatmap(mat, 
            display_numbers = F,
            color=colorRampPalette(c("white","black"))(10)[1:7],
            cluster_cols=F,cluster_rows = T,
            fontsize_row = 5,fontsize_col = 5)
save_plot(p,file.path(outdir_plots,"sup5.D.EMT_CAF_overlap.png"),width = 5,height = 11)


# EMT Epithelial
genes<-unique(unlist(c(top_ovary$PC2,top_breast$PC2,top_lung$PC2)))
mat<-cbind(genes %in% top_ovary$PC2,genes %in% top_breast$PC2,genes %in% top_lung$PC2)
rownames(mat)<-genes
colnames(mat)<-c("ovary","breast","lung")

row_sum<-rowSums(mat)
mat<-mat[order(row_sum,decreasing = T),]
mat[mat ==T]<-1
mat[mat ==F]<-0
p<-pheatmap(mat, 
            display_numbers = F,
            color=colorRampPalette(c("white","black"))(10)[1:7],
            cluster_cols=F,cluster_rows = T,
            fontsize_row = 5,fontsize_col = 5)
save_plot(p,file.path(outdir_plots,"sup5.D.EMT_Tumor_overlap.png"),width = 5,height = 11)




# sup5.e UMAP MAYA based with estrogen activity for 2 modes and corresponding violin plot with activity by cell type

mode="HALLMARK_ESTROGEN_RESPONSE_EARLY_mode1"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup5.E.UMAP_estrogen_mode1.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"sup5.E.Vln_estrogen_mode1_stat.png"),width = 14,height = 8)

mode="HALLMARK_ESTROGEN_RESPONSE_EARLY_mode2"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup5.E.UMAP_estrogen_mode2.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"sup5.E.Vln_estrogen_mode2_stat.png"),width = 14,height = 8)




# sup5.f UMAP MAYA based with coagulation activity for 3 modes and corresponding violin plot with activity by cell type


mode="HALLMARK_COAGULATION_mode1"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup5.F.UMAP_coagulation_mode1.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"sup5.F.Vln_coagulation_mode1.png"),width = 14,height = 8)

mode="HALLMARK_COAGULATION_mode2"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup5.F.UMAP_coagulation_mode2.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"sup5.F.Vln_coagulation_mode2.png"),width = 14,height = 8)

mode="HALLMARK_COAGULATION_mode3"
p<-plot_umap_gene(ovary_hallmark$umap,expr_mat = act_mat,gene = mode) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup5.F.UMAP_coagulation_mode3.png"),width = 7,height = 7)
p<-plot_violin_with_stat(act_mat,avg,mode=mode,meta_ovary,"Cell_type")
save_plot(p,file.path(outdir_plots,"sup5.F.Vln_coagulation_mode3.png"),width = 14,height = 8)


# sup5.g Heatmap metaprograms

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

# lines to adapt
p <- ggplot(data = nmf_intersect_melt_tumor, aes(x=Var1, y=Var2, fill=100*value/(100-value))) + 
    geom_tile() + 
    scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
    scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
    theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 6), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
    geom_vline(xintercept=c(424,504), linetype="longdash", size=0.6)+
    guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
save_plot(p,file.path(outdir_plots,"sup5.G.NMF_correlation_mat.png"),width = 15,height = 15)


# sup5.h Overlap MAYA/metaprograms

NMF_enrichment$pathway_prog<-paste0(NMF_enrichment$program, " - ", NMF_enrichment$pathway)
NMF_enrichment$pathway_prog<-factor(NMF_enrichment$pathway_prog,levels=rev(NMF_enrichment$pathway_prog))
p<- ggplot(NMF_enrichment,aes(x=pathway_prog,y=log10qval,fill=program))+
    geom_bar(stat="identity",width=0.8)+ 
    theme_classic()+
    scale_fill_manual(values=generate_palette(15)[10:15])+
    coord_flip()+
    theme(legend.position = "none",axis.title.y = element_blank())
save_plot(p,file= file.path(outdir_plots,"sup5.H.NMF_enrichment.png"),width = 20,height = 10)



# sup5.i Overlap MAYA/metaprograms

matrix_pct_tmp<-matrix_pct[rowMaxs(matrix_pct)>15,]

p<-pheatmap(matrix_pct_tmp,cluster_rows = T,cluster_cols = T,fontsize_row = 7,fontsize_col = 7,
            clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
            clustering_method = "ward.D2",color = magma(10,direction = -1),
            show_colnames = T,show_rownames = T)
save_plot(p,file.path(outdir_plots,"sup5.I.NMF_percentage_with_modes_filtered.png"),width = 13,height = 16)



###########



rm(list=ls())
gc()











