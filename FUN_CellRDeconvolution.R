function (Bulk, SingleCell, Bulk1, Alamat, MinCell, Mingene,
          Dimension, Alpha, maxi, mini)
{
  library(dplyr)
  library(Seurat)
  library(edgeR)
  library(glmnet)
  library(Matrix)
  pbmc.data <- SingleCell
  pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.genes = Mingene,
                             project = "CellR")
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data),
                     value = TRUE)
  percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes,
  ])/Matrix::colSums(pbmc@raw.data)
  pbmc <- AddMetaData(object = pbmc, metadata = percent.mito,
                      col.name = "percent.mito")
  pbmc <- AddMetaData(object = pbmc, metadata = Alamat, col.name = Cell)
  par(mfrow = c(1, 2))
  pbmc <- FilterCells(object = pbmc, subset.names = c("nGene",
                                                      "percent.mito"), low.thresholds = c(mini, -Inf), high.thresholds = c(maxi,
                                                                                                                           0.05))
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
                        scale.factor = 10000)
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean,
                            dispersion.function = LogVMR, do.plot = FALSE, x.low.cutoff = 0.0125,
                            x.high.cutoff = 3, y.cutoff = 0.5)
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI",
                                                       "percent.mito"))
  pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes,
                 do.print = FALSE, pcs.print = 1:5, genes.print = 5)
  pbmc <- SetAllIdent(object = pbmc, id = "Cell")
  Clusters <- pbmc@ident
  Levels <- levels(Clusters)
  Levels <- as.data.frame(Levels)
  Clusters <- as.data.frame(Clusters)
  Mahsool.Clusters <- Clusters
  cluster.markers <- matrix(0, 1, 6)
  cluster.markers <- as.data.frame(cluster.markers)
  colnames(cluster.markers) <- c("p_val", "avg_logFC", "pct.1",
                                 "pct.2", "p_val_adj", "cluster")
  cluster.markers <- cluster.markers[-c(1), ]
  for (i in 1:nrow(Levels)) {
    kk <- Levels[i, 1]
    kk <- as.data.frame(kk)
    kkk <- kk[1, 1]
    kkk <- as.character(kkk)
    print("Finding markers of the cluster", quote = FALSE)
    print(kkk, quote = FALSE)
    cluster1.markers <- FindMarkers(object = pbmc, ident.1 = Levels[i,
                                                                    1], min.pct = 0.25)
    a <- matrix(Levels[i, 1], nrow(cluster1.markers), 1)
    a <- as.data.frame(a)
    cluster1.markers$cluster <- a[, 1]
    cluster.markers <- rbind(cluster.markers, cluster1.markers)
  }
  sss <- which(cluster.markers$avg_logFC > 0)
  sss <- as.data.frame(sss)
  FinalMarkers <- cluster.markers[sss[, 1], ]
  wqwq <- which(rownames(FinalMarkers) %in% rownames(Bulk))
  wqwq <- as.data.frame(wqwq)
  FinalMarkers <- FinalMarkers[wqwq[, 1], ]
  cc <- pbmc@scale.data
  cc <- as.data.frame(cc)
  ssMM <- cc[rownames(FinalMarkers), ]
  ssMM <- as.data.frame(ssMM)
  rownames(ssMM) <- rownames(FinalMarkers)
  colnames(ssMM) <- colnames(cc)
  bulkTamiz <- cpm(Bulk, normalized.lib.sizes = TRUE, log = FALSE)
  bulkTamiz <- as.data.frame(bulkTamiz)
  bulkTamizFinal <- bulkTamiz[rownames(FinalMarkers), ]
  bulkTamizFinal <- as.data.frame(bulkTamizFinal)
  rownames(bulkTamizFinal) <- rownames(FinalMarkers)
  GTExNormalized <- Bulk1
  Mean <- rowMeans(Bulk1)
  Mean <- as.data.frame(Mean)
  Temp <- rep(Mean, nrow(GTExNormalized), ncol(GTExNormalized))
  Temp <- as.data.frame(Temp)
  Mian <- (GTExNormalized - Temp)
  Mian <- as.data.frame(Mian)
  rownames(Mian) <- rownames(GTExNormalized)
  colnames(Mian) <- colnames(GTExNormalized)
  STD <- matrix(0, nrow(GTExNormalized), 1)
  STD <- as.data.frame(STD)
  for (i in 1:nrow(GTExNormalized)) {
    STD[i, 1] <- sd(Mian[i, ], na.rm = TRUE)
  }
  Zero <- matrix(1e-05, nrow(STD), 1)
  Zero <- as.data.frame(Zero)
  STD <- Zero + STD
  Jam <- 0
  for (k in 1:nrow(STD)) {
    Jam <- (Jam + (1/STD[k, 1]))
  }
  Final <- matrix(0, nrow(GTExNormalized), 1)
  Final <- as.data.frame(Final)
  for (h in 1:nrow(STD)) {
    Final[h, 1] <- 1 + (1/(STD[h, 1] * Jam))
  }
  rownames(Final) <- rownames(GTExNormalized)
  ind <- matrix(0, nrow(ssMM), 1)
  for (i in 1:nrow(ssMM)) {
    if (is.na(ssMM[i, 1]) == TRUE) {
      ind[i, 1] <- 1
    }
  }
  ind <- as.data.frame(ind)
  www <- which(ind[, c(1)] == 0)
  www <- as.data.frame(www)
  FinalssMM <- ssMM[www[, 1], ]
  Finalbulk <- bulkTamizFinal[rownames(FinalssMM), ]
  Finalbulk <- as.data.frame(Finalbulk)
  rownames(Finalbulk) <- rownames(FinalssMM)
  size <- ncol(Finalbulk)
  Finalbulki <- matrix(0, nrow(Finalbulk), 1)
  Finalbulki <- as.data.frame(Finalbulki)
  rownames(Finalbulki) <- rownames(Finalbulk)
  GtexTemp <- Final[rownames(FinalssMM), ]
  GtexTemp <- rep(GtexTemp, 1, ncol(FinalssMM))
  FinalssMM <- as.matrix(FinalssMM)
  UniqueClusters <- unique(Clusters)
  row.names(UniqueClusters) <- UniqueClusters[, 1]
  OutputFinal <- matrix(0, nrow(UniqueClusters), ncol(Finalbulk))
  OutputFinal <- as.data.frame(OutputFinal)
  rownames(OutputFinal) <- row.names(UniqueClusters)
  OutputFreq <- OutputFinal
  for (Xha in 1:ncol(Finalbulk)) {
    Finalbulki <- as.data.frame(Finalbulki)
    Finalbulki <- Finalbulk[, Xha]
    Finalbulki <- as.matrix(Finalbulki)
    fit = glmnet(FinalssMM, Finalbulki, family = "gaussian",
                 alpha = Alpha, nlambda = 100)
    CEO <- coef(fit, s = 1)
    ttttt <- colnames(Finalbulki)
    ttttt <- as.data.frame(ttttt)
    Covariates <- matrix(0, nrow(Clusters), ncol(Finalbulki))
    Covariates <- as.data.frame(Covariates)
    rownames(Covariates) <- rownames(Clusters)
    colnames(Covariates) <- colnames(Finalbulki)
    tmp <- as.vector(CEO)
    tmp <- tmp[2:length(tmp)]
    Covariates[, 1] <- tmp
    UniqueClusters <- unique(Clusters)
    row.names(UniqueClusters) <- UniqueClusters[, 1]
    Output <- matrix(0, nrow(UniqueClusters), ncol(Covariates))
    Output <- as.data.frame(Output)
    row.names(Output) <- rownames(UniqueClusters)
    colnames(Output) <- colnames(Covariates)
    Frequency_Matrix <- Output
    for (i in 1:ncol(Covariates)) {
      Tempo <- Covariates[, i]
      Tempo <- as.data.frame(Tempo)
      rownames(Tempo) <- rownames(Covariates)
      Tempo <- cbind(Tempo, rownames(Covariates))
      TemInd <- which(Tempo[, 1] > 0)
      Tempo <- Tempo[c(TemInd), ]
      Tempo[, 2] <- Clusters[rownames(Tempo), 1]
      oo <- as.data.frame(table(Tempo[, 2]))
      rownames(oo) <- oo$Var1
      Percent <- oo$Freq/sum(oo$Freq) * 100
      oo$Var1 <- Percent
      Output[, i] <- oo[rownames(Output), 1]
      Frequency_Matrix[, i] <- oo[rownames(Output), c("Freq")]
    }
    OutputFinal[, Xha] <- Output
    OutputFreq[, Xha] <- Frequency_Matrix
  }
  colnames(OutputFinal) <- colnames(Bulk)
  colnames(OutputFreq) <- colnames(Bulk)
  FinalOutput <- list(Proportion = OutputFinal, Clusters = Clusters,
                      Frequency = OutputFreq)
  return(list(Proportion = OutputFinal, Markers = FinalMarkers,
              Clusters = Clusters, Frequency = OutputFreq))
}
