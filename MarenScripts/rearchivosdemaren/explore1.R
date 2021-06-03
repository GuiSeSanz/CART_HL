explore1 <- function(i, seurat_list, directoriodeseado){
  nm<-names(seurat_list)[[i]]
  print(nm)
  FIGURES_DIR=directoriodeseado
  print(file.path(FIGURES_DIR))
  
  seurat_work <- get (nm)
  seurat_work <- NormalizeData(object =  seurat_work,
                               normalization.method = "LogNormalize",scale.factor = 10000)
  seurat_work <- FindVariableFeatures(object = seurat_work,
                                      selection.method = "vst",nfeatures = 2000)
  top20 <- head(x = VariableFeatures(object = seurat_work),n =20)  # Identify the 20 most highly variable genes
  pdf(file.path(FIGURES_DIR,paste0("VariableTop20Features",nm,"_PRENORM.pdf")), width=11, height=8.5) 
  plot1 <- VariableFeaturePlot(object = seurat_work) # Plot variable features with labels
  print(LabelPoints(plot = plot1,points = top20, repel = TRUE))
  dev.off()
  
  plot1 <- LabelPoints(plot = plot1,points = top20, repel = TRUE)
  all_genes <- rownames(x = seurat_work)
  seurat_work <- ScaleData(object = seurat_work,features = all_genes) # Scale data
  seurat_work <- CellCycleScoring(seurat_work,g2m.features = g2m_genes,s.features = s_genes) # Perform cell cycle scoring 
  
  zz <- file(file.path(FIGURES_DIR,paste0("PCA_s",nm,"_PRENORM.txt")), open="wt") # Perform PCA and color by cell cycle phase
  sink(zz)
  sink(zz, type = "message")
  seurat_work <- RunPCA(seurat_work)
  sink(type="message")
  sink()
  closeAllConnections()
  
  pdf(file.path(FIGURES_DIR,paste0("PCA_bycellcycle_",nm,"_PRENORM.pdf")), width=11, height=8.5)
  print(DimPlot(seurat_work,reduction = "pca",group.by= "Phase", split.by = "Phase")) # Visualize the PCA, grouping by cell cycle phase
  dev.off()
  
  plot2<-DimPlot(seurat_work,reduction = "pca",group.by= "Phase")
  
  ggplot_lists <- Filter(function(x) is(x, "ggplot"), mget(ls()))
  print(paste0("Summary_",nm,"_PRENORM.pdf"))
  pdf(file.path(FIGURES_DIR,paste0("Summary_",nm,"_PRENORM.pdf")), width=11, height=8.5, onefile = T)
  #print(CombinePlots(ggplot_lists))
  print(plot1/plot2)
  dev.off()
  rm (plot1, plot2, ggplot_lists)
  
  ######ATENTTION I WONT REGRESS OUT THE SCORES ANY MORE, AFTER getting sucesfully substracting the "molt13" OK? 
  #vars_to_regress <- c("nUMI", "S.Score", "G2M.Score", "mitoRatio") # Define variables in metadata to regress
  vars_to_regress <- c("nUMI", "mitoRatio") # Define variables in metadata to regress
  # Regress out the uninteresting sources of variation in the data
  seurat_work <- ScaleData(object = seurat_work,vars.to.regress = vars_to_regress, verbose = FALSE) 
  seurat_work <- RunPCA(object = seurat_work)# Re-run the PCA
  
  pdf(file.path(FIGURES_DIR,paste0("PCA_bycellcycle_",nm,"_POSTSCALE.pdf")), width=11, height=8.5)
  print(DimPlot(object = seurat_work, reduction = "pca",group.by = "Phase", split.by = "Phase"))
  dev.off()
  
  pdf(file.path(FIGURES_DIR,paste0("PCA_Heatmap",nm,"_1.pdf")), width=11, height=8.5)
  print(DimHeatmap(seurat_work, dims = 1:15, cells = 500, balanced = TRUE)) # Explore heatmap of PCs
  dev.off()
  
  # Printing out the most variable genes driving PCs
  zz <- file(file.path(FIGURES_DIR,paste0("PCA_s",nm,"_1.txt")))
  sink(zz)
  sink(zz, type = "message")
  print(x = seurat_work[["pca"]], dims = 1:20, nfeatures = 10)
  sink(type="message")
  sink(NULL)
  closeAllConnections()
  
  ElbowPlot(object = seurat_work, ndims = 30) -> plot3 # Plot the elbow plot
  pct <- seurat_work[["pca"]]@stdev / sum(seurat_work[["pca"]]@stdev) * 100 # Determine percent of variation associated with each PC
  cumu <- cumsum(pct) # Calculate cumulative percents for each PC
  co1 <- which(cumu > 90 & pct < 5)[1] #PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  #print(co1)
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  #print(co2)# last point where change of % of variation is more than 0.1%.
  pcs <- min(co1, co2) #Calc Minimum of the two
  #print(pcs)
  
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct)) # Create a dataframe with values
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw() -> plot4 # Custom Elbow plot to visualize
  
  #print(x = seurat_work[["pca"]], dims = 1:25, nfeatures = 10) # Printing out the most variable genes driving PCs
  ggplot_lists <- Filter(function(x) is(x, "ggplot"), mget(ls()))
  print(paste0("Elbowplots",nm,".pdf"))
  pdf(file.path(FIGURES_DIR,paste0("Elbowplots",nm,".pdf")), width=8.5, height=11, onefile = T)
  #print(CombinePlots(ggplot_lists))
  print(plot3/plot4)
  dev.off()
  rm (plot3, plot4, ggplot_lists)
  assign(nm,seurat_work, envir = .GlobalEnv)
  #return()
}
#lapply(seq_along(LIST), explore1, seurat_list=LIST, directoriodeseado=DIR)

