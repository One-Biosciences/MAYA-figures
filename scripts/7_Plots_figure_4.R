### Plots Figure 4 ###

indir_files="./output_dir/"
outdir_plots="./output_plots/Figure4/"
dir.create(outdir_plots)

source("./scripts/utils.R")

### Load data  ###

# Kidney
umap_kidney_gene<-readRDS(file.path(indir_files,"umap_kidney_gene_based.rds"))
kidney_panglao<-readRDS(file.path(indir_files,"kidney_panglao.rds"))
meta_kidney<-readRDS(file.path(indir_files,"meta_kidney.rds"))
df_F1_kidney<-read.table(file=file.path(indir_files,"kidney_F1_scores.tsv"),sep="\t",header=T)

# CRC
umap_crc_gene<-readRDS(file.path(indir_files,"umap_crc_gene_based.rds"))
crc_panglao<-readRDS(file.path(indir_files,"crc_panglao.rds"))
meta_crc<-readRDS(file.path(indir_files,"meta_crc.rds"))
df_F1_crc<-read.table(file=file.path(indir_files,"CRC_F1_scores.tsv"),sep="\t",header=T)

# Ovary
umap_ovary_gene<-readRDS(file.path(indir_files,"umap_ovary_gene_based.rds"))
ovary_panglao<-readRDS(file.path(indir_files,"ovary_panglao.rds"))
meta_ovary<-readRDS(file.path(indir_files,"meta_ovary.rds"))

# Larynx
umap_larynx_gene<-readRDS(file.path(indir_files,"umap_larynx_gene_based.rds"))
larynx_panglao<-readRDS(file.path(indir_files,"larynx_panglao.rds"))
meta_larynx<-readRDS(file.path(indir_files,"meta_larynx.rds"))
final_df_patient_larynx<-read.table(file=paste0(indir_files,"larynx_shannon_results.tsv"),sep="\t",header=T)

# Pancreas
umap_panc8<-readRDS(file.path(indir_files,"umap_pancreas_gene_based.rds"))
panc8_panglao<-readRDS(file.path(indir_files,"pancreas_panglao.rds"))
meta_panc8<-readRDS(file.path(indir_files,"meta_panc8.rds"))
final_df_Technology<-read.table(file=paste0(indir_files,"panc8_shannon_results.tsv"),sep="\t",header=T)



###


# 4.a UMAP gene-based by author annotation kidney

p<-plot_umap_annot(umap_kidney_gene,type="Cell type",labels=meta_kidney[,"Cell_type"]) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"4.A.kidney_UMAP_author_annot_gene_based.png"),width = 7,height = 7)


# 4.b confusion matrix kidney

compare<-compare_clustering(meta_kidney$Cell_type,kidney_panglao$cell_annotation)

p<-pheatmap(compare,fontsize_row = 8,display_numbers = F,color=colorRampPalette(c("white","black"))(10)[1:9],main="MAYA prediction vs Author annotation",cluster_rows = F,cluster_cols = F)
save_plot(p,file.path(outdir_plots,"4.B.kidney_heatmap_MAYA_vs_author.png"),width = 15,height = 15)


# 4.c Boxplot F1 score kidney

{
    F1_plot<- df_F1_kidney[which(df_F1_kidney$metrics=="F1"),]
    
    citeCol <- generate_palette(5)
    names(citeCol) <- unique(meta_kidney$Cell_type)
    
    {
        symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", ""))
        
        p<-ggboxplot(
            as_tibble(F1_plot),
            x = "methods",
            y= "value",
            width = 0.5
        ) + geom_point(aes(fill = cell_type, x = methods, y = value), shape = 21, size = 4,position = position_jitterdodge(dodge.width = 0.4))+
            ggtitle("F1 score")+
            theme(
                aspect.ratio = 1,
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, face = "bold"),
                plot.title = element_text(hjust = 0.5)
            ) + grids(linetype = "dashed")+ scale_fill_manual(values = citeCol) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
            guides(color =F) +
            scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))+
            stat_compare_means(label = "p.signif", method = "wilcox.test",
                               ref.group = "MAYA",symnum.args=symnum.args,size=3)
    }
}
save_plot(p,file.path(outdir_plots,"4.C.kidney_F1_score.png"),width = 7,height = 7)


# 4.d UMAP gene-based by author annotation colon

p<-plot_umap_annot(umap_crc_gene,type="Cell type",labels=meta_crc[,"Cell_subtype"]) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"4.D.crc_UMAP_author_annot_gene_based.png"),width = 7,height = 7)


# 4.e confusion matrix colon

compare<-compare_clustering(meta_crc$Cell_subtype,crc_panglao$cell_annotation)
compare<-compare[,c("Enterocytes","Goblet cells","Pericytes","Smooth muscle cells","Dendritic cells","Macrophages","NK cells","T memory cells","B cells","Mast cells")]

p<-pheatmap(compare,fontsize_row = 8,display_numbers = F,color=colorRampPalette(c("white","black"))(10)[1:9],main="MAYA prediction vs Author annotation",cluster_rows = F,cluster_cols = F)
save_plot(p,file.path(outdir_plots,"4.E.CRC_heatmap_MAYA_vs_author.png"),width = 15,height = 15)


# 4.f Boxplot F1 score colon

{
    
    F1_plot<- df_F1_crc[which(df_F1_crc$metrics=="F1"),]
    
    citeCol <- generate_palette(10)
    names(citeCol) <- unique(meta_crc$Cell_subtype)
    
    {
        p<-ggboxplot(
            as_tibble(F1_plot),
            x = "methods",
            y= "value",
            width = 0.5
        ) + geom_point(aes(fill = cell_type, x = methods, y = value), shape = 21, size = 4,position = position_jitterdodge(dodge.width = 0.4))+
            ggtitle("F1 score")+
            theme(
                aspect.ratio = 1,
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, face = "bold"),
                plot.title = element_text(hjust = 0.5)
            ) + grids(linetype = "dashed")+ scale_fill_manual(values = citeCol) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
            guides(color =F) +
            scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))+
            stat_compare_means(label = "p.signif", method = "wilcox.test",
                               ref.group = "MAYA",symnum.args=symnum.args,size=3)
    }
}
save_plot(p,file.path(outdir_plots,"4.F.CRC_F1_score.png"),width = 7,height = 7)


# 4.g UMAP gene-based and MAYA based for larynx, by patient and cell type

# UMAP cell type gene based
p<-plot_umap_annot(umap_larynx_gene,type="Cell type",labels = meta_larynx$Cell_type,title = "Larynx - Author annotation - Gene based",colors = c(generate_palette(10)[c(3:6,9:10)]))
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"4.G.larynx_UMAP_author_annotation_no_legend.png"),width = 7,height = 7)

# UMAP patient gene based
p<-plot_umap_annot_random(umap_larynx_gene,type="Patient",labels = meta_larynx$patient_id,title = "Larynx - By patient - Gene based")
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"4.G.larynx_UMAP_patient_id__no_legend_random.png"),width = 7,height = 7)

# UMAP cell type MAYA based
p<-plot_umap_annot(larynx_panglao$umap,type="Cell type",labels = meta_larynx$Cell_type,title = "Larynx - Author annotation - Activity based",colors = c(generate_palette(10)[c(3:6,9:10)]))
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"4.G.larynx_UMAP_author_annotation_MAYA_no_legend.png"),width = 7,height = 7)

# UMAP patient MAYA based
p<-plot_umap_annot_random(larynx_panglao$umap,type="Patient",labels = meta_larynx$patient_id,title = "Larynx - By patient - Activity based")
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"4.G.larynx_UMAP_patient_MAYA_no_legend.png"),width = 7,height = 7)


# 4.h Boxplot SDI larynx


p<-ggboxplot(
    as_tibble(final_df_patient_larynx),
    x = "method",
    y= "shannon",
    width = 0.5
) + geom_point(aes(x = method, y = shannon,fill=method),shape = 21, size = 5,position = position_jitterdodge())+
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(hjust = 0.5)
    ) + grids(linetype = "dashed")+ggtitle("Shannon index by cluster for different methods")+
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.1))+
    scale_fill_manual(values=rep(wes_palette("IsleofDogs1")[6],times=4))+
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = "MAYA",symnum.args=symnum.args,size=3)

save_plot(p,file= file.path(outdir_plots,"4.H.larynx_shannon_index.png"),width = 15,height = 11)




# Supp

outdir_plots="./output_plots/Figure4_sup/"
dir.create(outdir_plots)


# sup4.a Boxplot recall accuracy kidney

{
    precision_plot<- df_F1_kidney[which(df_F1_kidney$metrics=="Precision"),]
    
    citeCol <- generate_palette(5)
    names(citeCol) <- unique(meta_kidney$Cell_type)
    
    {
        p<-ggboxplot(
            as_tibble(precision_plot),
            x = "methods",
            y= "value",
            width = 0.5
        ) + geom_point(aes(fill = cell_type, x = methods, y = value), shape = 21, size = 4,position = position_jitterdodge(dodge.width = 0.4))+
            ggtitle("Precision")+
            theme(
                aspect.ratio = 1,
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, face = "bold"),
                plot.title = element_text(hjust = 0.5)
            ) + grids(linetype = "dashed")+ scale_fill_manual(values = citeCol) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
            guides(color =F) +
            scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))+
            stat_compare_means(label = "p.signif", method = "wilcox.test",
                               ref.group = "MAYA",symnum.args=symnum.args,size=3)
    }
}
save_plot(p,file.path(outdir_plots,"sup4.A.kidney_precision_score.png"),width = 7,height = 7)

{
    recall_plot<- df_F1_kidney[which(df_F1_kidney$metrics=="Recall"),]
    
    citeCol <- generate_palette(5)
    names(citeCol) <- unique(meta_kidney$Cell_type)
    
    {
        p<-ggboxplot(
            as_tibble(recall_plot),
            x = "methods",
            y= "value",
            width = 0.5
        ) + geom_point(aes(fill = cell_type, x = methods, y = value), shape = 21, size = 4,position = position_jitterdodge(dodge.width = 0.4))+
            ggtitle("Recall")+
            theme(
                aspect.ratio = 1,
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, face = "bold"),
                plot.title = element_text(hjust = 0.5)
            ) + grids(linetype = "dashed")+ scale_fill_manual(values = citeCol) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
            guides(color =F) +
            scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))+
            stat_compare_means(label = "p.signif", method = "wilcox.test",
                               ref.group = "MAYA",symnum.args=symnum.args,size=3)
    }
}
save_plot(p,file.path(outdir_plots,"sup4.A.kidney_recall_score.png"),width = 7,height = 7)


# sup4.b Boxplot recall accuracy colon

{
    
    precision_plot<- df_F1_crc[which(df_F1_crc$metrics=="Precision"),]
    
    citeCol <- generate_palette(10)
    names(citeCol) <- unique(meta_crc$Cell_subtype)
    
    {
        p<-ggboxplot(
            as_tibble(precision_plot),
            x = "methods",
            y= "value",
            width = 0.5
        ) + geom_point(aes(fill = cell_type, x = methods, y = value), shape = 21, size = 4,position = position_jitterdodge(dodge.width = 0.4))+
            ggtitle("Precision")+
            theme(
                aspect.ratio = 1,
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, face = "bold"),
                plot.title = element_text(hjust = 0.5)
            ) + grids(linetype = "dashed")+ scale_fill_manual(values = citeCol) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
            guides(color =F) +
            scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))+
            stat_compare_means(label = "p.signif", method = "wilcox.test",
                               ref.group = "MAYA",symnum.args=symnum.args,size=3)
    }
}
save_plot(p,file.path(outdir_plots,"sup4.B.CRC_precision_score.png"),width = 7,height = 7)


{
    
    recall_plot<- df_F1_crc[which(df_F1_crc$metrics=="Recall"),]
    
    citeCol <- generate_palette(10)
    names(citeCol) <- unique(meta_crc$Cell_subtype)
    
    {
        p<-ggboxplot(
            as_tibble(recall_plot),
            x = "methods",
            y= "value",
            width = 0.5
        ) + geom_point(aes(fill = cell_type, x = methods, y = value), shape = 21, size = 4,position = position_jitterdodge(dodge.width = 0.4))+
            ggtitle("Recall")+
            theme(
                aspect.ratio = 1,
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, face = "bold"),
                plot.title = element_text(hjust = 0.5)
            ) + grids(linetype = "dashed")+ scale_fill_manual(values = citeCol) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
            guides(color =F) +
            scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "black"))+
            stat_compare_means(label = "p.signif", method = "wilcox.test",
                               ref.group = "MAYA",symnum.args=symnum.args,size=3)
    }
}
save_plot(p,file.path(outdir_plots,"sup4.B.CRC_recall_score.png"),width = 7,height = 7)


# sup4.c UMAP gene based for ovary

p<-plot_umap_annot(umap_ovary_gene,type="Cell type",labels = meta_ovary$Cell_type) + theme(legend.position = "none",plot.title=element_blank())
save_plot(p,file.path(outdir_plots,"sup4.C.Ovary_UMAP_author_annotation_gene_based.png"),width = 7,height = 7)


# sup4.d Confusion matrix ovary annotation

compare<-compare_clustering(meta_ovary$Cell_type,ovary_panglao$cell_annotation)
compare<-compare[c("EOC","CAF","Mesothelial","Endothelial","T-cells","NK","ILC","B-cells","Plasma-cells","pDC","Mast-cells","DC","Macrophages"),c("Epithelial cells","Fibroblasts","Endothelial cells","T cells","NK cells","B cells","Plasma cells","Dendritic cells","Macrophages")]

p<-pheatmap(compare,fontsize_row = 8,display_numbers = F,color=colorRampPalette(c("white","black"))(10)[1:9],main="MAYA prediction vs Author annotation",cluster_rows = F,cluster_cols = F)
save_plot(p,file.path(outdir_plots,"sup4.D.Ovary_heatmap_MAYA_vs_author.png"),width = 15,height = 15)


# sup4.e UMAP Gene-based and MAYA-based colored by cell type and technology

colors_tech<-generate_palette(5)
colors_cell_type<-generate_palette(20)[6:18]

# UMAP cell type gene based
p<-plot_umap_annot(umap_panc8,type="Cell type",labels = meta_panc8$celltype,title = "Pancreas - Author annotation - Gene based",colors = colors_cell_type)
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"sup4.E.Pancreas_UMAP_author_annotation_no_legend.png"),width = 7,height = 7)

# UMAP patient gene based
p<-plot_umap_annot_random(umap_panc8,type="Technology",labels = meta_panc8$tech,title = "Pancreas - By Technology - Gene based",colors = colors_tech)
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"sup4.E.Pancreas_UMAP_Technology_id__no_legend_random.png"),width = 7,height = 7)

# UMAP cell type based
p<-plot_umap_annot(test_panc8$umap,type="Cell type",labels = meta_panc8$celltype,title = "Pancreas - Author annotation - Activity based",colors = colors_cell_type)
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"sup4.E.Pancreas_UMAP_author_annotation_MAYA_no_legend.png"),width = 7,height = 7)

# UMAP Technology based
p<-plot_umap_annot_random(test_panc8$umap,type="Technology",labels = meta_panc8$tech,title = "Pancreas - By Technology - Activity based", colors = colors_tech)
save_plot_no_legend_no_title(p,file= file.path(outdir_plots,"sup4.E.Pancreas_UMAP_Technology_MAYA_no_legend.png"),width = 7,height = 7)


# sup4.Boxplot SDI pancreas

p<-ggboxplot(
    as_tibble(final_df_Technology),
    x = "method",
    y= "shannon",
    width = 0.5
) + geom_point(aes(x = method, y = shannon,fill=method),shape = 21, size = 5,position = position_jitterdodge())+
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(hjust = 0.5)
    ) + grids(linetype = "dashed")+ggtitle("Shannon index by cluster for different methods")+
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.1))+
    scale_fill_manual(values=rep(wes_palette("IsleofDogs1")[6],times=4))+
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = "MAYA",symnum.args=symnum.args,size=3)

save_plot(p,file= file.path(outdir_plots,"sup4.F.Pancreas_shannon_index.png"),width = 15,height = 11)



###########



rm(list=ls())
gc()





