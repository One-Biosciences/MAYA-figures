### Plots Figure 3 ###

indir_files="./output_dir/"
outdir_plots="./output_plots/Figure3/"
dir.create(outdir_plots)

source("./scripts/utils.R")

#################################



### Input Colon KEGG and Reactome

meta_crc<-readRDS(file.path(indir_files,"meta_crc.rds"))
logcpm_crc<-readRDS(file.path(indir_files,"logcpm_crc.rds"))
order_cells<-data.frame(cells=colnames(logcpm_crc),subtype=meta_crc$Cell_subtype)
order_cells <- order_cells %>%
    mutate(category =  factor(subtype, levels = unique(meta_crc$Cell_subtype))) %>%
    arrange(category)

crc_kegg<-readRDS(file.path(indir_files,"PCA_obj_KEGG_crc.rds"))
activity_mat_kegg<-build_activity_mat(crc_kegg$PCA_obj,scaled = T)
pagoda2_kegg<-as.matrix(read.table(file.path(indir_files,"crc_kegg_pagoda2.tsv"),header = T,check.names = F))
aucell_kegg<-as.matrix(read.table(file.path(indir_files,"crc_kegg_full_aucell.tsv"),header = T,check.names = F))

crc_reactome<-readRDS(file.path(indir_files,"PCA_obj_reactome_crc.rds"))
activity_mat_reactome<-build_activity_mat(crc_reactome$PCA_obj,scaled = T)
pagoda2_reactome<-as.matrix(read.table(file.path(indir_files,"crc_reactome_pagoda2.tsv"),header = T,check.names = F))
aucell_reactome<-as.matrix(read.table(file.path(indir_files,"crc_reactome_full_aucell.tsv"),header = T,check.names = F))



# 3.a top contributing genes kegg adhesion molecules

module="KEGG_CELL_ADHESION_MOLECULES_CAMS"
logcpm<-logcpm_crc[,order_cells$cells]
meta<-meta_crc[order_cells$cells,]
p<-plot_top_contrib(expr_mat=logcpm, PCA_object=crc_kegg$PCA_obj, module=module, n = 10, meta = meta,annot_name = "Cell_subtype",cluster_cols = F)
save_plot(p,file.path(outdir_plots,"3.A.crc_heatmap_kegg_CAMS_top_genes.png"),width = 15,height = 8)


# 3.b UMAP MAYA colored by annotation and by kegg adhesion molecules

p<-plot_umap_annot(crc_kegg$umap,type="Cell type",labels=meta_crc[,"Cell_subtype"]) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.B.crc_UMAP_kegg_author_annot_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_kegg$umap,expr_mat = activity_mat_kegg,gene = "KEGG_CELL_ADHESION_MOLECULES_CAMS_mode1") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.B.crc_UMAP_CAMS_mode1_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_kegg$umap,expr_mat = activity_mat_kegg,gene = "KEGG_CELL_ADHESION_MOLECULES_CAMS_mode2") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.B.crc_UMAP_CAMS_mode2_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_kegg$umap,expr_mat = activity_mat_kegg,gene = "KEGG_CELL_ADHESION_MOLECULES_CAMS_mode3") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.B.crc_UMAP_CAMS_mode3_no_legend.png"),width = 7,height = 7)


# 3.c top contributing genes reactome ion channels

module="REACTOME_ION_CHANNEL_TRANSPORT"
logcpm<-logcpm_crc[,order_cells$cells]
meta<-meta_crc[order_cells$cells,]
p<-plot_top_contrib(expr_mat=logcpm, PCA_object=crc_reactome$PCA_obj, module=module, n = 10, meta = meta,annot_name = "Cell_subtype",cluster_cols = F)
save_plot(p,file.path(outdir_plots,"3.C.crc_heatmap_reactome_ion_channel_transport_top_genes.png"),width = 15,height = 9)



# 3.d UMAP MAYA colored by annotation and by reactome ion channels

p<-plot_umap_annot(crc_reactome$umap,type="Cell type",labels=meta_crc[,"Cell_subtype"]) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.D.crc_UMAP_reactome_author_annot_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_reactome$umap,expr_mat = activity_mat_reactome,gene = "REACTOME_ION_CHANNEL_TRANSPORT_mode1") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.D.crc_UMAP_ION_mode1_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_reactome$umap,expr_mat = activity_mat_reactome,gene = "REACTOME_ION_CHANNEL_TRANSPORT_mode2") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.D.crc_UMAP_ION_mode2_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_reactome$umap,expr_mat = activity_mat_reactome,gene = "REACTOME_ION_CHANNEL_TRANSPORT_mode3") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.D.crc_UMAP_ION_mode3_no_legend.png"),width = 7,height = 7)

p<-plot_umap_gene(crc_reactome$umap,expr_mat = activity_mat_reactome,gene = "REACTOME_ION_CHANNEL_TRANSPORT_mode4") + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"3.D.crc_UMAP_ION_mode4_no_legend.png"),width = 7,height = 7)



### Supp

outdir_plots="./output_plots/Figure3_sup/"
dir.create(outdir_plots)

# sup3.a Activity matrix kegg colon

tmp<-activity_mat_kegg
rownames(tmp)<-sapply(rownames(tmp), function(x) strsplit(x,"KEGG_")[[1]][2])

p<-plot_heatmap_activity_mat(tmp,meta = meta_crc,annot_name = "Cell_subtype",cluster_rows = T,fontsize_row = 1.9)
save_plot(p,file.path(outdir_plots,"sup3.A.crc_heatmap_kegg_activity_matrix_scaled_subtype.png"),width = 15,height = 9)


# sup3.b Heatmap comparing MAYA modes to AUCell and Pagoda2 for kegg adhesion molecules

mat=matrix(c(crc_kegg$activity_matrix["KEGG_CELL_ADHESION_MOLECULES_CAMS_mode1",],
             crc_kegg$activity_matrix["KEGG_CELL_ADHESION_MOLECULES_CAMS_mode2",],
             crc_kegg$activity_matrix["KEGG_CELL_ADHESION_MOLECULES_CAMS_mode3",],
             aucell_kegg["KEGG_CELL_ADHESION_MOLECULES_CAMS",],
             pagoda2_kegg["KEGG_CELL_ADHESION_MOLECULES_CAMS",]),nrow=5,ncol=ncol(logcpm_crc),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2","MAYA_mode_3","AUCell","Pagoda2")
colnames(mat)<-colnames(activity_mat_kegg)

# order by subtype
mat<-mat[,order_cells$cells]

# Data frame with column annotations.
mat_col <- data.frame(Cell_type = meta_crc[order_cells$cells,]$Cell_subtype)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(10))
names(mat_colors$Cell_type) <- unique(meta_crc$Cell_subtype)

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
    main              = "CRC - KEGG_CELL_ADHESION_MOLECULES_CAMS",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8,
    gaps_col = cumsum(table(meta_crc$Cell_subtype)[unique(meta_crc$Cell_subtype)]),gaps_row = c(3,4)
)
save_plot(p,file.path(outdir_plots,"sup3.B.crc_heatmap_activity_kegg_cell_adh.png"),width = 15,height = 8)



# sup3.c Activity matrix reactome colon

tmp<-activity_mat_reactome
rownames(tmp)<-sapply(rownames(tmp), function(x) strsplit(x,"REACTOME_")[[1]][2])

p<-plot_heatmap_activity_mat(tmp,meta = meta_crc,annot_name = "Cell_subtype",cluster_rows = T,fontsize_row = 1.9)
save_plot(p,file.path(outdir_plots,"sup3.C.crc_heatmap_reactome_activity_matrix_scaled_subtype.png"),width = 20,height = 11)


# sup3.d Heatmap comparing MAYA modes to AUCell and Pagoda2 for reactome ion channels

mat=matrix(c(crc_reactome$activity_matrix["REACTOME_ION_CHANNEL_TRANSPORT_mode1",],
             crc_reactome$activity_matrix["REACTOME_ION_CHANNEL_TRANSPORT_mode2",],
             crc_reactome$activity_matrix["REACTOME_ION_CHANNEL_TRANSPORT_mode3",],
             crc_reactome$activity_matrix["REACTOME_ION_CHANNEL_TRANSPORT_mode4",],
             aucell_reactome["REACTOME_ION_CHANNEL_TRANSPORT",],
             pagoda2_reactome["REACTOME_ION_CHANNEL_TRANSPORT",]),nrow=6,ncol=ncol(logcpm_crc),byrow = T)
mat=scale_0_1(mat)
rownames(mat)<-c("MAYA_mode_1","MAYA_mode_2","MAYA_mode_3","MAYA_mode_4","AUCell","Pagoda2")
colnames(mat)<-colnames(logcpm_crc)

# order by subtype
mat<-mat[,order_cells$cells]

# Data frame with column annotations.
mat_col <- data.frame(Cell_type = meta_crc[order_cells$cells,]$Cell_subtype)
rownames(mat_col) <- colnames(mat)

mat_colors <- list(Cell_type = generate_palette(10))
names(mat_colors$Cell_type) <- unique(meta_crc$Cell_subtype)

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
    main              = "CRC - REACTOME_ION_CHANNEL_TRANSPORT",
    cluster_rows = F,
    cluster_cols = F,
    fontsize_row = 8,
    gaps_col = cumsum(table(meta_crc$Cell_subtype)[unique(meta_crc$Cell_subtype)]),
    gaps_row = c(4,5)
)
save_plot(p,file.path(outdir_plots,"sup3.D.crc_heatmap_activity_reactome_ion_channel.png"),width = 15,height = 8)


rm(list=ls())


