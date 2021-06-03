
explore2 <- function(i, seurat_list, seurat_dims_list, directoriodeseado, Resolution_Identities){
  nm<-names(seurat_list)[[i]]
  print(nm)
  FIGURES_DIR=directoriodeseado
  print(file.path(FIGURES_DIR))
  
  seurat_work <- get (nm)
  dimensions <- eval(parse(text=unname(seurat_dims_list[nm])))
  print(class(dimensions))
  
  Resolution_Identities<- Resolution_Identities
  for (i in seq_along(Resolution_Identities)) {
    print(Resolution_Identities[i])
    Idents(object = seurat_work) <- Resolution_Identities[i] #Assign identity of clusters
    set.seed(123)
    seurat_work <- RunTSNE(object = seurat_work)# Calculation of t-SNE
    plot5 <- DimPlot(object = seurat_work,label = TRUE,reduction = "tsne")  # Plotting t-SNE
    plot5 <- plot5 + ggtitle("t-SNE")
    set.seed(123)
    seurat_work <- RunUMAP(seurat_work, reduction = "pca", dims = dimensions) # Calculation of UMAP
    plot6 <- DimPlot(seurat_work,reduction = "umap",label = TRUE,label.size = 6) # Plot the UMAP
    plot6 <- plot6 + ggtitle("UMAP")
    ggplot_lists <- Filter(function(x) is(x, "ggplot"), mget(ls()))
    print(paste0("Explore_Tsne_UMAP",nm,"_",Resolution_Identities[i],".pdf"))
    pdf(file.path(FIGURES_DIR,paste0("Explore_Tsne_UMAP",nm,"_",Resolution_Identities[i],".pdf")),
        width=11, height=8.5, onefile = T)
    #print(CombinePlots(ggplot_lists))
    print(plot5|plot6)
    rm (plot5, plot6, ggplot_lists)
    dev.off() #https://github.com/satijalab/seurat/issues/2697
    # Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
    n_cells <- FetchData(seurat_work, vars = c("ident")) %>% dplyr::count(ident) %>% spread(ident, n)
    
    # View table
    # View(n_cells)  #####I could print this as well and capture it in a unique .txt
    group_by <- c("Phase") # Establishing groups to color plots by
    
    # Getting coordinates for cells to use for UMAP and associated grouping variable information
    class_umap_data <- FetchData(seurat_work, vars = c("ident", "UMAP_1", "UMAP_2", group_by))
    umap_label <- FetchData(seurat_work, vars = c("ident", "UMAP_1", "UMAP_2"))  %>% group_by(ident) %>% 
      summarise(x=mean(UMAP_1), y=mean(UMAP_2)) # Adding cluster label to center of cluster on UMAP
    
    # Getting coordinates for cells to use for PCA and associated grouping variable information
    class_pca_data <- FetchData(seurat_work, vars = c("ident", "PC_1", "PC_2", group_by))
    pca_label <- FetchData(seurat_work, vars = c("ident", "PC_1", "PC_2"))  %>% mutate(ident = seurat_work@active.ident) %>%
      group_by(ident) %>% summarise(x=mean(PC_1), y=mean(PC_2)) # Adding cluster label to center of cluster on PCA
    # Function to plot UMAP and PCA as grids
    pdf(file.path(FIGURES_DIR,paste0("UMAP_PCAs",nm,"_",Resolution_Identities[i],".pdf")),
        width=11, height=8.5, onefile = T)
    plot_grid(
      ggplot(class_umap_data, aes(UMAP_1, UMAP_2)) +
        geom_point(aes_string(color = group_by), alpha = 0.7) +
        scale_color_brewer(palette = "Set2")  +
        geom_text(data=umap_label, aes(label=ident, x, y)),
      ggplot(class_pca_data, aes(PC_1, PC_2)) +
        geom_point(aes_string(color = group_by), alpha = 0.7) +
        scale_color_brewer(palette = "Set2")  +
        geom_text(data=pca_label, 
                  aes(label=ident, x, y)),
      nrow = 1, 
      align = "v") -> plot7
    print(plot7)
    dev.off()
    
    # Determine metrics to plot present in seurat_control@meta.data
    metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
    # Extract the UMAP coordinates for each cell and include information about the metrics to plot
    qc_data <- FetchData(seurat_work, vars = c(metrics, "ident", "UMAP_1", "UMAP_2"))
    
    # Plot a UMAP plot for each metric
    pdf(file.path(FIGURES_DIR,paste0("UMAP_for_each_metric",nm,"_",Resolution_Identities[i],".pdf")),
        width=11, height=8.5, onefile = T)
    map(metrics, function(qc){
      ggplot(qc_data,
             aes(UMAP_1, UMAP_2)) +
        geom_point(aes_string(color=qc), 
                   alpha = 0.7) +
        scale_color_gradient(guide = FALSE, 
                             low = "grey90", 
                             high = "blue")  +
        geom_text(data=umap_label, 
                  aes(label=ident, x, y)) +
        ggtitle(qc)
    }) %>%
      plot_grid(plotlist = .) -> plot8
    print(plot8)
    dev.off()
    
    # Defining the information in the seurat object of interest
    columns <- c(paste0("PC_", dimensions),"ident","UMAP_1", "UMAP_2")
    pc_data <- FetchData(seurat_work, vars = columns) # Extracting this data from the seurat object
    
    # Plotting a UMAP plot for each of the PCs
    pdf(file.path(FIGURES_DIR,paste0("PCAs_for_each_metric",nm,"_",Resolution_Identities[i],".pdf")),
        width=11, height=8.5, onefile = T)
    map(paste0("PC_", dimensions), function(pc){
      ggplot(pc_data, 
             aes(UMAP_1, UMAP_2)) +
        geom_point(aes_string(color=pc), 
                   alpha = 0.7) +
        scale_color_gradient(guide = FALSE, 
                             low = "grey90", 
                             high = "blue")  +
        geom_text(data=umap_label, 
                  aes(label=ident, x, y)) +
        ggtitle(pc)
    }) %>% 
      plot_grid(plotlist = .) -> plot9
    print(plot9)
    dev.off()
    rm(plot7, plot8, plot9)
    
    pdf(file.path(FIGURES_DIR,paste0("UMAPS_For_KnownMarkers",nm,"_",Resolution_Identities[i],".pdf")),
        width=11, height=8.5, onefile = T)
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("CD79A", "CD19","CD14", "CAR-pCCL-BCMA","CD3D", "PTPRC"), label = T))
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("CD4", "CD8A", "EOMES", "TBX21", "FOXP3"), label = T))
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("SELL","CCR7","CXCR3", "CD27", "CD28", "CD127"), label = T))

    print(FeaturePlot(seurat_work, reduction = "umap", features = c("IL7R","KLF2","TCF7"), label = T))
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("GZMA", "GZMB", "GZMK", "GZMH", "GNLY", "NKG7"), label = T))
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("HLA-DRA", "IL2RA", "IL3RA", "CD69", "CD40LG"), label = T))

    print(FeaturePlot(seurat_work, reduction = "umap", features = c("E2F4", "E2F7", "CDK1"), label = T))
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("TIGIT", "LAG3","PDCD1","CTLA4"), label = T))
    print(FeaturePlot(seurat_work, reduction = "umap", features = c("TOX", "TOX2", "IRF8", "IRF4", "NR4A2", "ENTPD1"), label = T))

    dev.off()
    
    pdf(file.path(FIGURES_DIR,paste0("ViolinPlot_KnownMarkers",nm,"_",Resolution_Identities[i],".pdf")),
        width=11, height=8.5, onefile = T)
    print(VlnPlot(seurat_work, features = c("CD79A", "CD19","CD14", "CAR-pCCL-BCMA","CD3D", "PTPRC"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("CD4", "CD8A", "EOMES", "TBX21", "FOXP3"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("SELL","CCR7","CXCR3", "CD27", "CD28", "CD127"), pt.size = -1))

    print(VlnPlot(seurat_work, features = c("CXCR3", "CD27", "CD28"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("FOXP3","GNLY", "NKG7"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("CD44","CREM","CD69"), pt.size = -1))

    print(VlnPlot(seurat_work, features = c("IL7R","KLF2","TCF7"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("GZMA", "GZMB", "GZMK", "GZMH", "GNLY", "NKG7"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("HLA-DRA", "IL2RA", "IL3RA", "CD69", "CD40LG"), pt.size = -1))

    print(VlnPlot(seurat_work, features = c("E2F4", "E2F7", "CDK1"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("TIGIT", "LAG3","PDCD1","CTLA4"), pt.size = -1))
    print(VlnPlot(seurat_work, features = c("TOX", "TOX2", "IRF8", "IRF4", "NR4A2", "ENTPD1"), pt.size = -1))
   
    dev.off()
    #####
    #####
  }
  
}
