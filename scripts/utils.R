
### Libraries ###

library(msigdbr)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(Seurat)
library(SingleCellExperiment)
library(MAYA)
library(ggpubr)
library(SCINA)
library(SeuratData)
library(AUCell)
library(pagoda2)
library(CellID)
library(harmony)
library(scales)
library(RColorBrewer)
library(viridis)
library(fgsea)


### variables ###

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", ""))


### Metrics ###

compare_clustering<-function(test_annot,ref_annot){
    ref_clusters<-unique(ref_annot)
    test_clusters<-unique(test_annot)
    mat<-matrix(0,ncol=length(ref_clusters),nrow=length(test_clusters))
    colnames(mat)<-unique(ref_clusters)
    rownames(mat)<-unique(test_clusters)
    
    for(i in 1:length(test_clusters)){
        for(j in 1:length(ref_clusters)){
            tmp_test<-which(test_annot==test_clusters[i])
            tmp_ref<-which(ref_annot==ref_clusters[j])
            mat[i,j]<-round(length(intersect(tmp_test,tmp_ref))/length(tmp_test),2)
        }
    }
    
    return(mat)
}

average_by_cluster<-function(mat,annotation){
    
    # get all possible clusters
    clusters<-unique(annotation)
    
    # compute row means for each subsetted matrix by cell type and store it as a list
    means<-lapply(clusters,function(x){
        if(length(which(annotation==x))>1){
            tmp<-mat[,which(annotation==x)]
            rowMeans(tmp)
        }
        else{
            mat[,which(annotation==x)]
        }
        
    })
    # store result as a matrix
    names(means)<-clusters
    out<-as.matrix(dplyr::bind_cols(means))
    rownames(out)<-rownames(mat)
    
    return(out)
}

compute_SDI<-function(annotation,cluster_annot){
    # numerize annotation
    annotation_binary<-rep(NA,times=length(annotation))
    i=1
    for(annot in unique(annotation)){
        annotation_binary[which(annotation==annot)]<-i
        i<-i+1
    }
    # compute Shannon Diversity Index
    clusters<-c()
    shannon<-c()
    for(cluster in unique(cluster_annot)){
        clusters<-c(clusters,cluster)
        prop<-prop.table(table(as.numeric(annotation_binary)[which(cluster_annot==cluster)]))
        prop<-prop*log(prop)
        shannon<-c(shannon,(-1)*sum(prop)/log(length(unique(as.numeric(annotation_binary)))))
    }
    df<-data.frame(cluster=clusters,shannon=shannon)
    df[order(df$cluster),]
    return(df)
}


compute_F1_score<-function(predictions="prediction",reference="Cell_type",meta = meta_kidney,mapping = map_cite_kidney,data ="kidney"){
    
    TruthDF <-
        tibble(
            Author = as.vector(meta[,reference]),
            Predictions = as.vector(meta[,predictions]),
            map_label = as.vector(sapply(mapping[as.vector(meta[,reference])], function(x)
                paste0(x, collapse = ", ")))
        )
    #TruthDF <- TruthDF %>% dplyr::filter(map_label != "removed")
    TruthDF <-
        TruthDF %>% dplyr::mutate(positive = mapply(
            x = .$map_label,
            y = .$Predictions,
            FUN = function(x, y) {
                x %like% y
            }
        ))
    TruthDF <-
        TruthDF %>% dplyr::mutate(final_map = ifelse(positive, Author, Predictions))
    
    PredDF <- sapply(unique(TruthDF$Author), function(x, truth) {
        TruthDFpos <- TruthDF %>% dplyr::filter(Author == x)
        TruthDFneg <- TruthDF %>% dplyr::filter(Author != x)
        TP <- sum(TruthDFpos$Author == TruthDFpos$final_map)
        FN <- sum(TruthDFpos$Author != TruthDFpos$final_map)
        FP <- sum(TruthDFneg$final_map %in% mapping[[x]])
        Recall <- TP / (TP + FN)
        Precision <- TP / (FP + TP)
        F1 <- 2 * (Recall * Precision) / (Recall + Precision)
        return(c(
            Precision = Precision,
            Recall = Recall,
            F1 = F1
        ))
    },
    truth = TruthDF) %>% as.data.frame() %>%  rownames_to_column(var = "metrics") %>%  gather("cell_type", "value",-1) %>%  dplyr::mutate(methods = predictions) %>%  dplyr::mutate(data = data) %>%  dplyr::mutate(value = ifelse(is.na(value), 0, value))
    
    return(PredDF)
}



### Plots ###

plot_umap_annot_random<-function(umap,type="Cell type",labels,title=NULL,colors=NULL){
    
    #### Build plot ####
    df <- data.frame(x = umap$layout[,1],
                     y = umap$layout[,2],
                     Alias = labels)
    if(is.null(colors)){
        palette=generate_palette(length(unique(labels)))
        p<-ggplot(df[sample(1:nrow(df),size = nrow(df),replace = F),], aes(x, y, colour = Alias)) +
            geom_point(shape=21,colour="white",aes(fill=Alias),size=2,stroke=0.5)+
            xlab("UMAP_1")+ylab("UMAP_2")+
            ggtitle(title)+theme_classic()+
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values=palette,breaks=unique(labels)) +
            labs(fill = type)
    }
    else{
        p<-ggplot(df[sample(1:nrow(df),size = nrow(df),replace = F),], aes(x, y, colour = Alias)) +
            geom_point(shape=21,colour="white",aes(fill=Alias),size=2,stroke=0.5)+
            xlab("UMAP_1")+ylab("UMAP_2")+
            ggtitle(title)+theme_classic()+
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values=colors,breaks=unique(labels)) +
            labs(fill = type)
    }
    print(p)
}



plot_top_contrib<-function (expr_mat, PCA_object, module, n = 10, meta = NULL,annot_name = NULL, cluster_cols = T, fontsize = 7, colors_annot = NULL,scale=T)
{
    stopifnot(is.list(PCA_object), is.character(module), is.numeric(n),
              n > 1)
    if (!(module %in% names(PCA_object))) {
        stop(paste0("No informative activity available for ",
                    module))
    }
    nb_comp <- length(PCA_object[[module]]$expl_var)
    top_genes <- c()
    contrib <- t(PCA_object[[module]]$gene_contrib)
    if (nb_comp == 1) {
        contrib <- contrib[order(contrib, decreasing = T), ]
        names <- names(contrib)
        top_genes <- names[1:n]
    }
    else {
        comp <- sapply(colnames(contrib), function(x) strsplit(x,
                                                               "PC")[[1]][2])
        list_top_genes <- lapply(comp, function(x) {
            contrib <- contrib[order(contrib[, paste0("PC",
                                                      x)], decreasing = T), ]
            names <- names(contrib[, paste0("PC", x)])
            names[1:n]
        })
        top_genes <- unlist(list_top_genes)
    }
    if (!is.null(annot_name)) {
        mat_col <- data.frame(as.factor(meta[, annot_name]))
        colnames(mat_col) <- annot_name
        rownames(mat_col) <- colnames(expr_mat)
        mat_colors <- list()
        if (!is.null(colors_annot)) {
            if (length(unique(meta[, annot_name])) == length(colors_annot[[annot_name]])) {
                mat_colors[[annot_name]] <- colors_annot[[annot_name]]
                names(mat_colors[[annot_name]]) <- unique(meta[,
                                                               annot_name])
            }
            else {
                stop(paste0("Number or colors provided in colors_annot different from number of levels of ",
                            annot_name))
            }
        }
        else {
            palette = generate_palette(length(unique(meta[,
                                                          annot_name])))
            names(palette) = unique(meta[, annot_name])
            mat_colors <- list(palette)
            names(mat_colors) <- annot_name
        }
    }
    if (nb_comp > 1) {
        gaps_row = seq(from = n, to = n * nb_comp - n, by = n)
    }
    else {
        gaps_row = n
    }
    
    if(scale==T){
        mat = MAYA::scale_0_1(Matrix::t(scale(Matrix::t(expr_mat[top_genes, ]),
                                              center = T, scale = T)))
    }
    else{
        mat = expr_mat[top_genes, ]
    }
    
    if (!is.null(annot_name)) {
        p <- pheatmap(mat = mat, color = inferno(100, direction = 1)[20:100],
                      border_color = NA, show_colnames = FALSE, show_rownames = TRUE,
                      annotation_col = mat_col, annotation_colors = mat_colors,
                      drop_levels = TRUE, fontsize = fontsize, main = paste0("Top contributor genes - ",
                                                                             module), cluster_rows = F, cluster_cols = cluster_cols,
                      gaps_row = gaps_row, clustering_distance_cols = "euclidean",
                      clustering_method = "ward.D2")
    }
    else {
        p <- pheatmap(mat = mat, color = inferno(100, direction = 1)[20:100],
                      border_color = NA, show_colnames = FALSE, show_rownames = TRUE,
                      drop_levels = TRUE, fontsize = fontsize, main = paste0("Top contributor genes - ",
                                                                             module), cluster_rows = F, cluster_cols = cluster_cols,
                      gaps_row = gaps_row, clustering_distance_cols = "euclidean",
                      clustering_method = "ward.D2")
    }
}


plot_violin_with_stat<-function(act_mat,avg,mode,meta,annot_name,colors=NULL){
    symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", ""))
    df<-data.frame(cells=colnames(act_mat),score=act_mat[mode,])
    df$Cell_type<-meta[df$cells,annot_name]
    df$Cell_type<-factor(df$Cell_type, levels =unique(df$Cell_type))
    ref=colnames(avg)[which.max(avg[mode,])]
    n=length(unique(df$Cell_type))
    
    if(is.null(colors)){
        colors=generate_palette(n)
    }
    
    ggplot(df,aes(x=Cell_type,y=score,fill=Cell_type,colour=Cell_type))+geom_violin()+
        theme_classic()+scale_fill_manual(values = colors,breaks=unique(df$Cell_type))+
        scale_color_manual(values = colors,breaks=unique(df$Cell_type))+
        theme(legend.position = "none",axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+
        ylim(c(0,1))+stat_compare_means(label = "p.signif", method = "wilcox.test",
                                        ref.group = ref,symnum.args=symnum.args,size=3)
}

### Saving plots ###

save_plot<-function(p,file,res=300,width=7,height=7){
    png(file,res = res,width = width,height = height,units="cm")
    print(p)
    dev.off()
    paste0("Plot saved to: ", file)
}


save_plot_no_legend_no_title<-function(p,file,res=300,width=7,height=7){
    png(file,res = res,width = width,height = height,units="cm")
    print(p + theme(legend.position = "none",plot.title =element_blank()))
    dev.off()
    paste0("Plot saved to: ", file)
}



### Extract info from MAYA object ###


get_top_contrib<-function(PCA_object, module, n = 10){
    nb_comp <- length(PCA_object[[module]]$expl_var)
    top_genes <- list()
    contrib <- t(PCA_object[[module]]$gene_contrib)
    if (nb_comp == 1) {
        PC=names(PCA_object[[module]]$expl_var)
        contrib <- contrib[order(contrib, decreasing = T), ]
        names <- names(contrib)
        top_genes[[PC]] <- names[1:n]
    }
    else {
        comp <- sapply(colnames(contrib), function(x) strsplit(x, 
                                                               "PC")[[1]][2])
        list_top_genes <- lapply(comp, function(x) {
            contrib <- contrib[order(contrib[, paste0("PC", 
                                                      x)], decreasing = T), ]
            names <- names(contrib[, paste0("PC", x)])
            names[1:n]
        })
        top_genes <- list_top_genes
    }
    return(top_genes)
}

get_diff_contrib<-function(PCA_object, module){
    nb_comp <- length(PCA_object[[module]]$expl_var)
    max_contrib <- list()
    contrib <- t(PCA_object[[module]]$gene_contrib)
    if (nb_comp == 1) {
        contrib<-contrib[order(contrib,decreasing = T)]
        PC=names(PCA_object[[module]]$expl_var)
        max_contrib[[PC]] <- contrib[1]-contrib[2]
    }
    else {
        comp <- sapply(colnames(contrib), function(x) strsplit(x, 
                                                               "PC")[[1]][2])
        list_top_genes <- lapply(comp, function(x) {
            contrib <- contrib[order(contrib[, paste0("PC", 
                                                      x)], decreasing = T), ]
            tmp<-contrib[, paste0("PC", x)]
            tmp[1]-tmp[2]
        })
        max_contrib <- list_top_genes
    }
    return(max_contrib)
}


df_contrib_diff_ratio<-function(PCA_obj){
    comp<-c()
    diff_contrib<-c()
    max_contrib<-c()
    for(module in names(PCA_obj)){
        comp<-c(comp,names(PCA_obj[[module]][["expl_var"]]))
        diff_contrib<-c(diff_contrib,unlist(get_diff_contrib(PCA_obj,module=module)))
        max_contrib<-c(max_contrib,unlist(get_max_contrib(PCA_obj,module=module)))
    }
    df<-data.frame(comp,diff_contrib,max_contrib,ratio=diff_contrib/max_contrib)
    return(df)
}

get_nb_contrib<-function(PCA_object, module, thr = 0.05){
    nb_comp <- length(PCA_object[[module]]$expl_var)
    top_genes <- list()
    contrib <- t(PCA_object[[module]]$gene_contrib)
    if (nb_comp == 1) {
        PC=names(PCA_object[[module]]$expl_var)
        contrib <- contrib[order(contrib, decreasing = T), ]
        names <- names(contrib)
        top_genes[[PC]] <- names[which(contrib>=thr)]
    }
    else {
        comp <- sapply(colnames(contrib), function(x) strsplit(x, 
                                                               "PC")[[1]][2])
        list_top_genes <- lapply(comp, function(x) {
            contrib <- contrib[order(contrib[, paste0("PC", 
                                                      x)], decreasing = T), ]
            names <- names(contrib[, paste0("PC", x)])
            names[which(contrib[, paste0("PC", x)]>=thr)]
        })
        top_genes <- list_top_genes
    }
    return(top_genes)
}





### Pathway enrichment with hypergeometric test ###

load_MSIGdb <- function(ref, GeneSetClasses){
    if ((!ref %in% c("hg38", "mm10")) ) 
        stop("ChromSCape::gene_set_enrichment_analysis_scExp - ",
             "Reference genome (ref) must be ",
             "'hg38' or 'mm10' if gene sets not specified.")
    stopifnot(is.character(GeneSetClasses))
    message(
        paste0(
            "ChromSCape::gene_set_enrichment_analysis_scExp - Loading ",
            ref,
            " MSigDB gene sets."
        )
    )
    columns = c("gs_name", "gs_cat", "gene_symbol")
    if (ref == "hg38")
        GeneSetsDf = msigdbr::msigdbr("Homo sapiens")[, columns]
    if (ref == "mm10")
        GeneSetsDf = msigdbr::msigdbr("Mus musculus")[, columns]
    colnames(GeneSetsDf) = c("Gene.Set", "Class", "Genes")
    system.time({
        GeneSetsDf <- GeneSetsDf %>% dplyr::group_by(
            .data$Gene.Set, .data$Class) %>%
            dplyr::summarise("Genes" = paste(.data$Genes,
                                             collapse = ","))
    })
    corres = data.frame(
        long_name = c("c1_positional", "c2_curated", "c3_motif", 
                      "c4_computational", "c5_GO", "c6_oncogenic",
                      "c7_immunologic", "hallmark"), short_name = c(
                          paste0("C", seq_len(7)), "H"))
    GeneSetsDf$Class = corres$long_name[
        match(GeneSetsDf$Class, corres$short_name)]
    GeneSetsDf = GeneSetsDf[which(GeneSetsDf$Class %in% GeneSetClasses),]
    GeneSets = lapply(GeneSetsDf$Gene.Set, function(x) {
        unlist(strsplit(
            GeneSetsDf$Genes[which(GeneSetsDf$Gene.Set == x)], split = ","))})
    names(GeneSets) = GeneSetsDf$Gene.Set
    return(GeneSets)
}


enrichment_markers_provided_markers<-function(ref="hg38",GeneSetClasses=c("c2_curated","hallmark"),top_markers,qval=0.05){
    if(!exists("gene_sets_msigdb")){
        gene_sets_msigdb <- load_MSIGdb(ref,GeneSetClasses)
    }
    
    gene_sets=gene_sets_msigdb
    sep = ";"
    possibleIds <- unique(unlist(gene_sets))
    
    mylist <- unique(top_markers)
    gene.sets <- lapply(gene_sets, unique)
    nids <- length(possibleIds)
    gene.sets <- lapply(gene.sets, function(x) intersect(x, possibleIds))
    nref <- as.numeric(lapply(gene.sets, length))
    gene.sets <- gene.sets[nref > 0]
    n <- length(mylist)
    fun <- function(x) {
        y <- intersect(x, mylist)
        nx <- length(x)
        ny <- length(y)
        pval <- stats::phyper(ny - 1, nx, nids - nx, n, lower.tail = FALSE)
        c(nx, ny, pval, paste(y, collapse = sep))
    }
    tmp <- as.data.frame(t(as.matrix(vapply(gene.sets, fun, FUN.VALUE = c(
        "Nb_of_genes" = 0, "Nb_of_deregulated_genes" = 0,
        "p-value" = 0, "Deregulated_genes" = ""
    )))))
    rownames(tmp) <- names(gene.sets)
    for (i in seq_len(3)) {
        tmp[, i] <- as.numeric(
            as.character(tmp[, i])
        )
    }
    tmp <- data.frame(
        tmp[, seq_len(3)], p.adjust(tmp[, 3], method = "BH"), tmp[, 4]
    )
    names(tmp) <- c(
        "Nb_of_genes", "Nb_of_deregulated_genes",
        "p-value", "q-value", "Deregulated_genes"
    )
    tmp
    tmp<-tmp[which(tmp$`q-value`<qval),]
    tmp<-tmp[order(tmp$`q-value`),]
    
    return(tmp)
}



### Run Seurat CCA ###

run_seuratCCA<-function(seurat,annot_name){
    seurat.list <- SplitObject(seurat, split.by = annot_name)
    
    # perform standard preprocessing on each object
    for (i in 1:length(seurat.list)) {
        seurat.list[[i]] <- NormalizeData(seurat.list[[i]], verbose = FALSE)
        seurat.list[[i]] <- FindVariableFeatures(
            seurat.list[[i]], selection.method = "vst",
            nfeatures = 2000, verbose = FALSE
        )
    }
    
    # find anchors
    anchors <- FindIntegrationAnchors(object.list = seurat.list)
    
    # integrate data
    integrated <- IntegrateData(anchorset = anchors)
    
    
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated)
    integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))
    integrated <- RunUMAP(integrated, dims = 1:10)
    
    return(integrated)
}


### generate panglao list by organ if specified or all

generate_panglao_list<-function(organs=NULL){
    path<-system.file("extdata", "PanglaoDB_markers_27_Mar_2020.tsv", package = "MAYA")
    panglao<-read.table(path,header=T,sep="\t",quote='')
    panglao<-panglao[which(panglao$species=="Mm Hs" | panglao$species=="Hs"),]
    
    # basic types
    panglao_basic<-panglao[which(panglao$organ %in% c("Connective tissue","Smooth muscle","Immune system","Vasculature","Blood","Epithelium","Skeletal muscle")),]
    basic_panglao<-list()
    for(cell_type in unique(panglao_basic$cell.type)){
        if(!(cell_type %in% c("Endothelial cells (blood brain barrier)","Endothelial cells (aorta)","Gamma delta T cells"))){
            basic_panglao[[cell_type]]<-panglao_basic[which(panglao_basic$cell.type==cell_type),"official.gene.symbol"]
        }
        
    }
    
    ### add lists for other organs if required
    if(!is.null(organs)){
        if(!is.null(intersect(unique(panglao$organ),organs))){
            panglao_specif<-panglao[which(panglao$organ %in% organs),]
            
            for(cell_type in unique(panglao_specif$cell.type)){
                basic_panglao[[cell_type]]<-panglao_specif[which(panglao_specif$cell.type==cell_type),"official.gene.symbol"]
            }
            
        }
        if(organs=="all"){
            for(cell_type in unique(panglao$cell.type)){
                basic_panglao[[cell_type]]<-panglao[which(panglao$cell.type==cell_type),"official.gene.symbol"]
            }
            
        }
        
    }
    return(basic_panglao)
}


### GSEA ###

GSEA <- function(gene_list,GO_file,pval_thr=0.01){
    set.seed(54321)
    fgseaRes <- fgsea(pathways = GO_file, 
                      stats    = gene_list,
                      eps      = 0.0,
                      minSize  = 15,
                      maxSize  = 500)
    fgseaRes <- as.data.frame(fgseaRes) %>%  dplyr::filter(padj < !!pval_thr) 
    fgseaRes$logqval<-(-log10(fgseaRes$padj))
    fgseaRes <- fgseaRes %>% arrange(desc(NES))
    fgseaRes$up_down <- ifelse(fgseaRes$NES >= 0, "Up", "Down")
    return(fgseaRes)
}



### NMF ###

robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
    
    # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
    intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
    intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
    nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
    names(nmf_sel) <- names(nmf_programs)
    
    # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
    # ii) the minimum overlap with programs from another cell line
    nmf_sel_unlist <- do.call(cbind, nmf_sel)
    inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
    
    final_filter <- NULL 
    for(i in names(nmf_sel)) {
        a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
        b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
        if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 10
        if(length(b) > 1) {
            c <- names(b[1]) 
            for(y in 2:length(b)) {
                if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
            }
            final_filter <- c(final_filter, c)
        } else {
            final_filter <- c(final_filter, names(b))
        }
    }
    return(final_filter)                                                      
}


### MAYA without informativity filters ###

run_activity_analysis_no_filter_info<-function(expr_mat,modules_list,norm=FALSE,nb_comp_max=5,min_cells_pct=0.05,nCores=1,min_module_size=10,max_contrib=0.5){
    
    message("Computing gene sets activity")
    
    #### Check parameters ####
    stopifnot(is.list(modules_list),
              is.logical(norm),is.numeric(nb_comp_max),
              is.numeric(min_cells_pct),is.numeric(nCores),
              is.numeric(min_module_size),is.numeric(max_contrib))
    stopifnot(is(expr_mat, 'sparseMatrix') | is(expr_mat, 'matrix'))
    
    expr_mat<-as(expr_mat,"dgCMatrix")
    
    #### Normalize data if specified ####
    if(norm){
        expr_mat<-logcpmNormalization(expr_mat)
    }
    # remove genes expressed in less than 10 cells
    expr_mat<-expr_mat[Matrix::rowSums(expr_mat!=0)>=10,]
    
    ### get minimum number of cells ###
    min_cells<-round(min_cells_pct*ncol(expr_mat))
    
    #### Run PCA for each gene set ####
    compute_gene_set_activity <- function(x) {
        common_genes<-intersect(modules_list[[x]],rownames(expr_mat))
        if(length(common_genes)>=min_module_size){
            # perform PCA, compute explained var and modify PCA sign to favor activation.
            pca<-stats::prcomp(scale(Matrix::t(expr_mat[common_genes,]),center=T,scale=T), scale = FALSE,retx = TRUE,rank. = 20)
            n=ncol(pca$x)
            ExpVar <- apply(pca$x[, 1:n], 2, var)/sum(apply(scale(Matrix::t(expr_mat[common_genes,]),center=T,scale=T), 2, var))
            # orientate PCs to favor activation (change orientation if genes with negative contributions have higher weights in absolute value than positive one)
            for(j in 1:min(min_module_size,nb_comp_max)){
                if(sum(pca$rotation[,j])<0){
                    pca$x[,j]<-(-pca$x[,j])
                    pca$rotation[,j]<-(-pca$rotation[,j])
                }
            }
            # scale projection for threshold detection
            projection<-scale_0_1(t(pca$x))
            # assess informativity of successive PCs and stop when conditions are unmet.
            comp<-c()
            list_thr<-c()
            nb_selected_cells<-c()
            for(i in 1:min(min_module_size,nb_comp_max)){
                thr<-.activity_assignmentThreshold(activity = projection[i,])
                comp<-c(comp,i)
                list_thr<-c(list_thr,thr)
                nb_selected_cells<-c(nb_selected_cells,length(which(projection[i,]>=thr)))
            }
            if(length(comp)!=0){
                nb_comp<-min(length(comp),nb_comp_max,min_module_size)
                # keep only PCs that pass threshold and keep scaled score
                if(length(comp)==1){
                    activity_scores<-as.matrix(t(projection[comp,])) # cells as columns
                    rownames(activity_scores)<-paste0("mode",comp)
                }else{
                    activity_scores<-projection[comp,]
                    rownames(activity_scores)<-paste0("mode",comp)
                }
                activity_scores_raw<-t(pca$x[,comp])
                rownames(activity_scores_raw)<-paste0("mode",comp)
                gene_contrib<-t(pca$rotation[,comp]) # genes as columns
                expl_var<-ExpVar[comp]
                list(activity_scores=activity_scores,activity_scores_raw=activity_scores_raw,gene_contrib=gene_contrib,expl_var=expl_var,list_thr=list_thr,nb_selected_cells=nb_selected_cells)
            }
        }
    }
    
    
    output_object<-parallel::mclapply(names(modules_list), FUN = compute_gene_set_activity,mc.cores = nCores)
    
    
    #### Generate final output ####
    names(output_object)<-names(modules_list)
    #remove modules with no informative scores from the output object.
    output_object<-output_object[lapply(output_object,length)!=0]
    
    message("Found at least one informative activation mode for ",length(output_object)," gene sets")
    
    return(output_object)
}


