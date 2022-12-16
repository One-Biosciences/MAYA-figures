### Plots Figure 1 ###

indir_files="./output_dir/"
dir.create("./output_plots/")
outdir_plots="./output_plots/Figure1_sup/"
dir.create(outdir_plots)

source("./scripts/utils.R")

#################################



### boxplot variance explained by PC all pathways ###

no_info_kidney_hallmark<-readRDS(file=paste0(indir_files,"PCA_obj_no_info_hallmark_kidney.rds"))
no_info_colon_kegg<-readRDS(file=paste0(indir_files,"PCA_obj_no_info_KEGG_crc.rds"))
no_info_colon_reactome<-readRDS(file=paste0(indir_files,"PCA_obj_no_info_reactome_crc.rds"))


### hallmark

# boxplot variance explained by PC all pathways

pathway_list<-c()
var_explained<-c()
for(pathway in names(no_info_kidney_hallmark)){
    pathway_list<-c(pathway_list, rep(pathway,5))
    var_explained<-c(var_explained,no_info_kidney_hallmark[[pathway]]$expl_var)
}
mode<-rep(1:5,times=length(no_info_kidney_hallmark))

df<-data.frame(pathway=pathway_list,mode=as.character(mode),var_explained=var_explained)

p<-df %>%
    ggplot(aes(x=mode,y=var_explained, fill=mode)) +
    geom_boxplot(aes(group=mode),outlier.shape = NA) +
    geom_line(aes(group=pathway),color="lightgrey",position = position_dodge(0.2)) +
    geom_point(aes(fill=mode,group=pathway), position = position_dodge(0.2)) +
    theme_classic()+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values=generate_palette(10)[6:10])+
    ggtitle("Kidney - HALLMARK")+
    xlab("Mode")+
    ylab("Variance explained")

save_plot(p,file= file.path(outdir_plots,"sup1.A.boxplot_kidney_hallmark_var_explained_by_mode.png"),width = 15,height = 11)



# kegg


pathway_list<-c()
var_explained<-c()
for(pathway in names(no_info_colon_kegg)){
    pathway_list<-c(pathway_list, rep(pathway,5))
    var_explained<-c(var_explained,no_info_colon_kegg[[pathway]]$expl_var)
}
mode<-rep(1:5,times=length(no_info_colon_kegg))

df<-data.frame(pathway=pathway_list,mode=as.character(mode),var_explained=var_explained)

p<-df %>%
    ggplot(aes(x=mode,y=var_explained, fill=mode)) +
    geom_boxplot(aes(group=mode),outlier.shape = NA) +
    geom_line(aes(group=pathway),color="lightgrey",position = position_dodge(0.2),size=0.1) +
    geom_point(aes(fill=mode,group=pathway), position = position_dodge(0.2)) +
    theme_classic()+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values=generate_palette(10)[6:10])+
    ggtitle("Colon - KEGG")+
    xlab("Mode")+
    ylab("Variance explained")

save_plot(p,file= file.path(outdir_plots,"sup1.A.boxplot_colon_kegg_var_explained_by_mode_all.png"),width = 15,height = 11)


# reactome


pathway_list<-c()
var_explained<-c()
for(pathway in names(no_info_colon_reactome)){
    pathway_list<-c(pathway_list, rep(pathway,5))
    var_explained<-c(var_explained,no_info_colon_reactome[[pathway]]$expl_var)
}
mode<-rep(1:5,times=length(no_info_colon_reactome))

df<-data.frame(pathway=pathway_list,mode=as.character(mode),var_explained=var_explained)

p<-df %>%
    ggplot(aes(x=mode,y=var_explained, fill=mode)) +
    geom_boxplot(aes(group=mode),outlier.shape = NA) +
    geom_line(aes(group=pathway),color="lightgrey",position = position_dodge(0.2),size=0.1) +
    geom_point(aes(fill=mode,group=pathway), position = position_dodge(0.2)) +
    theme_classic()+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values=generate_palette(10)[6:10])+
    ggtitle("Colon - REACTOME")+
    xlab("Mode")+
    ylab("Variance explained")

save_plot(p,file= file.path(outdir_plots,"sup1.A.boxplot_colon_reactome_var_explained_by_mode_all.png"),width = 15,height = 11)





### SS2 vs 10X ###


# data 
pbmcsca_seurat<-readRDS(file.path(indir_files,"pbmcsca_seurat.rds"))
SS2_hallmark<-readRDS(file.path(indir_files,"pbmcsca_SS2_hallmark.rds"))
tenX_hallmark<-readRDS(file.path(indir_files,"pbmcsca_tenX_hallmark.rds"))


# VennDiagrams

pathways<-list(SS2=names(SS2_hallmark$PCA_obj),
               `10X`=names(tenX_hallmark$PCA_obj))

png(file.path(outdir_plots,"sup1.H.venn_diag_shared_pathways_SS2_10X.png"),res=300,height=5,width=7,units = "cm")
grid.newpage();
venn.plot <- draw.pairwise.venn(area1 = length(pathways$SS2),area2 = length(pathways$`10X`),cross.area = length(intersect(pathways$SS2,pathways$`10X`)),category = c("SS2","10X"),fill=c("#0073C2FF", "#EFC000FF"),alpha=0.2)
grid.draw(venn.plot)
dev.off()

fisher_test_from_list(pathways,50)
# 0.001240679

# shared modes for pathways in common: to see how many modes we found for common pathways with each modality (removing different pathways from diff numbers)
pathways_in_common<-intersect(names(SS2_hallmark$PCA_obj),names(tenX_hallmark$PCA_obj))

modes<-list(SS2=rownames(build_activity_mat(SS2_hallmark$PCA_obj[pathways_in_common])),
            `10X`=rownames(build_activity_mat(tenX_hallmark$PCA_obj[pathways_in_common])))

png(file.path(outdir_plots,"sup1.H.venn_diag_shared_modes_common_pathways_SS2_10X.png"),res=300,height=5,width=7,units = "cm")
grid.newpage();
venn.plot <- draw.pairwise.venn(area1 = length(modes$SS2),area2 = length(modes$`10X`),cross.area = length(intersect(modes$SS2,modes$`10X`)),category = c("SS2","10X"),fill=c("#0073C2FF", "#EFC000FF"),alpha=0.2)
grid.draw(venn.plot)
dev.off()

fisher_test_from_list(modes,length(pathways_in_common)*5)
# 1.837471e-23

# for the common modes: compute the number of overlapping genes between top10 contributors

top_contrib_SS2<-list()
top_contrib_10X<-list()
for(pathway in pathways_in_common){
    top_contrib_SS2[[pathway]]<-get_top_contrib(SS2_hallmark$PCA_obj,pathway,n=10)
    top_contrib_10X[[pathway]]<-get_top_contrib(tenX_hallmark$PCA_obj,pathway,n=10)
}

for(pathway in pathways_in_common){
    PC=intersect(names(top_contrib_SS2[[pathway]]),names(top_contrib_10X[[pathway]]))
    top_contrib_SS2[[pathway]]<-top_contrib_SS2[[pathway]][PC]
    top_contrib_10X[[pathway]]<-top_contrib_10X[[pathway]][PC]
}

size_overlap<-c()
for(pathway in pathways_in_common){
    n=length(top_contrib_SS2[[pathway]])
    if(n>0){
        for(i in 1:n){
            size_overlap<-c(size_overlap,length(intersect(top_contrib_SS2[[pathway]][[i]],top_contrib_10X[[pathway]][[i]])))
        }
    }
    
}
df<-data.frame(overlap=size_overlap,method=rep("",times=length(size_overlap)))

p<-gghistogram(
    as_tibble(df),
    x = "overlap",
    bins=11,fill="method",add="mean"
) + 
    theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(hjust = 0.5)
    )+ggtitle("Overlap between top10 contributing genes")+
    scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(0,10))+
    scale_fill_manual(values=rep(wes_palette("IsleofDogs1")[6],times=1))

save_plot(p,file= file.path(outdir_plots,"sup1.H.hist_number_overlaping_genes_top10.png"),width = 15,height = 11)


# Distributions examples 

# HALLMARK_TNFA_SIGNALING_VIA_NFKB
module="HALLMARK_TNFA_SIGNALING_VIA_NFKB"
thr=SS2_hallmark[["PCA_obj"]][[module]][["list_thr"]][1]

png(file.path(outdir_plots,"sup1.H.distrib_TNFA_SS2_10X.png"),res=300,height=7,width=15,units = "cm")

par(mfrow=c(1,2))
hist(SS2_hallmark[["PCA_obj"]][[module]][["activity_scores"]],breaks=100,main="Smart-Seq2",xlab="Scaled activity")
abline(v=thr,col="red", lwd=3, lty=2)

thr=tenX_hallmark[["PCA_obj"]][[module]][["list_thr"]][1]
hist(tenX_hallmark[["PCA_obj"]][[module]][["activity_scores"]][1,],breaks=100,main="10X (v3)",xlab="Scaled activity")
abline(v=thr,col="red", lwd=3, lty=2)
dev.off()

# HALLMARK_ALLOGRAFT_REJECTION_mode1
module="HALLMARK_ALLOGRAFT_REJECTION"
thr=SS2_hallmark[["PCA_obj"]][[module]][["list_thr"]][1]

png(file.path(outdir_plots,"sup1.H.distrib_ALLO_SS2_10X.png"),res=300,height=7,width=15,units = "cm")

par(mfrow=c(1,2))
hist(SS2_hallmark[["PCA_obj"]][[module]][["activity_scores"]],breaks=100,main="Smart-Seq2",xlab="Scaled activity")
abline(v=thr,col="red", lwd=3, lty=2)

thr=tenX_hallmark[["PCA_obj"]][[module]][["list_thr"]][1]
hist(tenX_hallmark[["PCA_obj"]][[module]][["activity_scores"]][1,],breaks=100,main="10X (v3)",xlab="Scaled activity")
abline(v=thr,col="red", lwd=3, lty=2)
dev.off()

# HALLMARK_INFLAMMATORY_RESPONSE_mode1
module="HALLMARK_INFLAMMATORY_RESPONSE"
thr=SS2_hallmark[["PCA_obj"]][[module]][["list_thr"]][1]

png(file.path(outdir_plots,"sup1.H.distrib_INFL_SS2_10X.png"),res=300,height=7,width=15,units = "cm")

par(mfrow=c(1,2))
hist(SS2_hallmark[["PCA_obj"]][[module]][["activity_scores"]],breaks=100,main="Smart-Seq2",xlab="Scaled activity")
abline(v=thr,col="red", lwd=3, lty=2)

thr=tenX_hallmark[["PCA_obj"]][[module]][["list_thr"]][1]
hist(tenX_hallmark[["PCA_obj"]][[module]][["activity_scores"]][1,],breaks=100,main="10X (v3)",xlab="Scaled activity")
abline(v=thr,col="red", lwd=3, lty=2)

dev.off()









