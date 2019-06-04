#
# Copyright (C) 2017 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

MLClusteringDensityBased <- function(jaspResults, dataset, options, ...) {
  
  .ss <- function(x) sum(scale(x, scale = FALSE)^2)
  .disttovar <- function(x) {
    mean(x**2)/2
  }
  .tss <- function(x) {
    n <- nrow(as.matrix(x))
    .disttovar(x)*(n-1)
  }
  

  # read variables ##
  dataset             <- .densityClusteringReadData(dataset, options)
  
  # error handling & code variable names in base64
  .densityClusteringErrorHandling(dataset, options)
  ready               <- length(options[["predictors"]][options[["predictors"]] != ""] > 0)
  
  # Run the analysis and save the results
  res                 <- .densityClustering(dataset, options, jaspResults, ready)
  
  # create the evaluation table
  .densityClusteringSummaryTable(res, options, jaspResults, ready)
  
  # create the cluster information table
  .densityClusteringInformationTable(options, res, jaspResults, ready)
  
  # Create the cluster plot
  .densityBasedClusterPlot(dataset, options, res, jaspResults, ready)
  
  # Create the k-distance plot
  .kdistPlot(options, res, jaspResults, ready)
}

.densityClusteringReadData <- function(dataset, options){
  predictors <- unlist(options[['predictors']])
  predictors <- predictors[predictors != ""]
  if (is.null(dataset)) {
    dataset <- .readDataSetToEnd(columns.as.numeric = predictors, exclude.na.listwise = predictors)
  }
  if(options[["scaleEqualSD"]]){
    dataset <- scale(dataset)
  }
  return(dataset)
}

.densityClusteringErrorHandling <- function(dataset, options){
  predictors <- unlist(options$predictors)
  if(length(predictors[predictors != '']) > 0){
    for(i in 1:length(predictors)){
      errors <- .hasErrors(dataset, perform, type = c('infinity', 'observations'),
                           all.target = predictors[i],
                           observations.amount = "< 2",
                           exitAnalysisIfErrors = TRUE)
    }
  }
}

.densityClustering <- function(dataset, options, jaspResults, ready){
  
  if(!is.null(jaspResults[["res"]]$object)) return(jaspResults[["res"]]$object)
  
  # set the seed so that every time the same set is chosen (to prevent random results) ##
  if(options[["seedBox"]])
    set.seed(options[["seed"]])
  
  if(ready){
    if(options[["modelOpt"]] == "validationManual"){
      res <- .densityClusteringManual(dataset, options)
    } else {
      res <- .densityClusteringOptimized(dataset, options)
    }
  } else {
    res <- list()
  }
  jaspResults[["res"]] <- createJaspState(res)
  jaspResults[["res"]]$dependOn(c("predictors", "eps", "minPts", "modelOpt", "seed", 
                                         "seedBox", "scaleEqualSD"))
  
  return(jaspResults[["res"]]$object)
}

.densityClusteringManual <- function(dataset, options){
  
  densityfit <- dbscan::dbscan(as.data.frame(dataset[, .v(options[["predictors"]])]),
                        eps = options[['eps']],
                        minPts = options[['minPts']])
  
  df <- as.data.frame(dataset[, .v(options[["predictors"]])])
  cs = fpc::cluster.stats(dist(df), densityfit$cluster)
  if(options[["distance"]] == "Correlated densities") {
    df <- as.dist(1-cor(dataset[, .v(options[["predictors"]])],
                          method="pearson"))
    cs <- fpc::cluster.stats(df, densityfit$cluster)
  } else {
    df = df
    cs = cs
  }

  res <- list()
  res[['Predictions']] <- data.frame(
    'Observation' = 1:nrow(dataset),
    'Cluster' = densityfit$cluster
  )
  res[["pred.values"]] <- densityfit$cluster
  res[["clusters"]] <- length(table(densityfit$cluster))-1
  res[['eps']] <- densityfit$eps
  res[["N"]] <- nrow(dataset)
  res[['size']] <- table(densityfit$cluster)
  res[["noisePoints"]] <- length(densityfit$cluster[densityfit$cluster == 0])
  res[["MinPts"]] <- densityfit$minPts
  res[["Silh_score"]] <- cs$avg.silwidth
  res[['WSS']] <- cs$within.cluster.ss
  res[['TSS']] <- .tss(df)
  res[['BSS']] <- res[['TSS']]-res[['WSS']]
  D = cs$within.cluster.ss
  data <- as.data.frame(dataset[, .v(options[["predictors"]])])
  m = ncol(kmeans(data, res[["clusters"]])$centers)
  n = length(densityfit$cluster)
  k = length(table(densityfit$cluster))-1
  res[['AIC']] <- D + 2*m*k
  res[['BIC']] <- D + log(n)*m*k
  res[['silh_scores']] <- cs$clus.avg.silwidths
  withinss <- 0
  for (i in 0:length(table(densityfit$cluster))) {
    withinss[i] <- .ss(data[densityfit$cluster == i])
  }
  res[['WSS']] <- withinss
  return(res)
}

.densityClusteringOptimized <- function(dataset, options){
  
  df <- as.data.frame(dataset[, .v(options[["predictors"]])])
  cs = fpc::cluster.stats(dist(df), densityfit$cluster)
  if(options[["distance"]] == "Correlated densities") {
    df <- as.dist(1-cor(dataset[, .v(options[["predictors"]])],
                        method="pearson"))
    cs <- fpc::cluster.stats(df, densityfit$cluster)
  } else {
    df = df
    cs = cs
  }
  
  minPts <- dim(df)[2]*2
  optics_eps <- dbscan::optics(df, minPts = minPts) 
  eps <- optics_eps$eps
  densityfit <- dbscan::dbscan(df, eps = eps, minPts = minPts)
  
  res <- list()
  res[['Predictions']] <- data.frame(
    'Observation' = 1:nrow(dataset),
    'Cluster' = densityfit$cluster
  )
  res[["pred.values"]] <- densityfit$cluster
  res[["clusters"]] <- length(table(densityfit$cluster))-1
  res[['eps']] <- densityfit$eps
  res[["N"]] <- nrow(dataset)
  res[['size']] <- table(densityfit$cluster)
  res[["noisePoints"]] <- length(densityfit$cluster[densityfit$cluster == 0])
  res[["MinPts"]] <- densityfit$minPts
  res[["Silh_score"]] <- cs$avg.silwidth
  res[['WSS']] <- cs$within.cluster.ss
  res[['TSS']] <- .tss(df)
  res[['BSS']] <- res[['TSS']]-res[['WSS']]
  D = cs$within.cluster.ss
  data <- as.data.frame(dataset[, .v(options[["predictors"]])])
  m = ncol(kmeans(data, res[["clusters"]])$centers)
  n = length(densityfit$cluster)
  k = length(table(densityfit$cluster))-1
  res[['AIC']] <- D + 2*m*k
  res[['BIC']] <- D + log(n)*m*k
  res[['silh_scores']] <- cs$clus.avg.silwidths
  withinss <- 0
  for (i in 0:length(table(densityfit$cluster))) {
     withinss[i] <- .ss(data[densityfit$cluster == i])
  }
  res[['WSS']] <- withinss
  return(res)
}
  
.densityClusteringSummaryTable <- function(res, options, jaspResults, ready){
  
  if(!is.null(jaspResults[["evaluationTable"]])) return() #The options for this table didn't change so we don't need to rebuild it
  
  evaluationTable                       <- createJaspTable("Density based Clustering Model Summary")
  jaspResults[["evaluationTable"]]      <- evaluationTable
  jaspResults[["evaluationTable"]]$position <- 1
  evaluationTable$dependOn(c("predictors", "eps", "minPts", "modelOpt", "seed", "scaleEqualSD"))
  
  evaluationTable$addColumnInfo(name = 'clusters', title = 'Clusters', type = 'integer')
  evaluationTable$addColumnInfo(name = 'eps', title = 'Eps', type = 'integer')
  evaluationTable$addColumnInfo(name = 'minPts', title = 'MinPts', type = 'integer')
  evaluationTable$addColumnInfo(name = 'noisePoints', title = 'Noise Points', type = 'integer')
  evaluationTable$addColumnInfo(name = 'measure', title = 'R\u00B2', type = 'number', format = 'dp:2')
  evaluationTable$addColumnInfo(name = 'aic', title = 'AIC', type = 'number', format = 'dp:1')
  evaluationTable$addColumnInfo(name = 'bic', title = 'BIC', type = 'number', format = 'dp:1')
  evaluationTable$addColumnInfo(name = 'Silh', title = 'Silhouette value', type = 'number', format = 'dp:1')
  evaluationTable$addColumnInfo(name = 'n', title = 'N', type = 'number', format = 'dp:1')
  
  if(!ready)
    return()
  
  row <- data.frame(clusters = res[["clusters"]], eps = res[["eps"]], minPts = res[["MinPts"]],
                    noisePOints = res[["noisePoints"]], n = round(res[["N"]], 0), measure = res[['BSS']]/res[['TSS']],
                    aic = res[['AIC']], bic = res[['BIC']], Silh = res[['Silh_score']])
  }

.densityClusteringInformationTable <- function(options, res, jaspResults, ready){
  
  if(!is.null(jaspResults[["clusterInfoTable"]])) return() #The options for this table didn't change so we don't need to rebuild it
  
  if (options[['tableClusterInformation']]){
    
    clusterInfoTable                        <- createJaspTable("Cluster Information")
    jaspResults[["clusterInfoTable"]]       <- clusterInfoTable
    clusterInfoTable$dependOn(c("tableClusterInformation","predictors", "eps", "minPts", "modelOpt",
                                "scaleEqualSD"))
    clusterInfoTable$position               <- 2
    clusterInfoTable$transpose              <- TRUE
    
    clusterInfoTable$addColumnInfo(name = 'cluster', title = 'Cluster', type = 'integer')
    clusterInfoTable$addColumnInfo(name = 'size', title = 'Size', type = 'integer')
    clusterInfoTable$addColumnInfo(name = 'withinss_table', title = 'Within Sum of Squares', type = 'number', format = 'dp:2')
    clusterInfoTable$addColumnInfo(name = 'silh_scores', title = 'Silhouette scores', type = 'integer')

    if(!ready)
      return()
    
    size <- res[["size"]]
    cluster <- 0:res[["clusters"]]
    withinss_table <- res[["WSS"]]
    silh_scores <- res[['silh_scores']]
    row <- data.frame(cluster = cluster, size = size, withinss_table = withinss_table,
                      silh_scores = silh_scores)

  }
}

.densityBasedClusterPlot <- function(dataset, options, res, jaspResults, ready){
  if(!ready && options[['plot2dCluster']]){
    p <- createJaspPlot(plot = NULL, title= "Cluster Plot", width = 400, height = 300)
    p$setError("Plotting not possible: No analysis has been run.")
    return()
  } else if(ready && options[['plot2dCluster']]){
    if(is.null(jaspResults[["plot2dCluster"]])){
      
      if(options[["seedBox"]])
        set.seed(options[["seed"]])
      
      unique.rows <- which(!duplicated(dataset[, .v(options[["predictors"]])]))
      data <- dataset[unique.rows, .v(options[["predictors"]])]
      tsne_out <- Rtsne::Rtsne(as.matrix(data))
      
      densityfit <- dbscan::dbscan(data, eps = res[["eps"]], minPts = res[['MinPts']])
      pred.values <- densityfit$cluster
      
      tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = pred.values)
      p <- ggplot2::ggplot(tsne_plot) + ggplot2::geom_point(ggplot2::aes(x = x, y = y, fill = factor(col)), size = 4, stroke = 1, shape = 21, color = "black") + ggplot2::xlab(NULL) + ggplot2::ylab(NULL)
      p <- p + ggplot2::scale_fill_manual(values = colorspace::qualitative_hcl(n = res[["clusters"]]))
      p <- JASPgraphs::themeJasp(p)
      p <- p + ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())
      jaspResults[["plot2dCluster"]] 		<- createJaspPlot(plot = p, title= "T-sne Cluster Plot", width = 400, height = 300)
      jaspResults[["plot2dCluster"]]		$dependOn(options =c("predictors", "eps", "minPts", "modelOpt", "seed",
                                                           "scaleEqualSD", "ready", "plot2dCluster"))
      jaspResults[["plot2dCluster"]] 		$position <- 3
    }
  }
}

.kdistPlot <- function(options, res, jaspResults, ready){
  if(!ready && options[['k-distplot']] && options[["modelOpt"]] != "validationManual"){
    p <- createJaspPlot(plot = NULL, title= "K-distance plot", width = 400, height = 300)
    p$setError("Plotting not possible: No analysis has been run.")
    return()
  } else if(ready && options[['k-distplot']] && options[["modelOpt"]] != "validationManual"){
    if(is.null(jaspResults[["optimPlot"]])){
      
      data <- dataset[unique.rows, .v(options[["predictors"]])]
      
      jaspResults[["optimPlot"]] 		 <- createJaspPlot(plot = JASPgraphs::themeJasp(kNNdist(data, k = res[["MinPts"]])),
                                                      title= "K-distance plot", height = 300, width = 400)
      jaspResults[["optimPlot"]]		   $dependOn(options =c("predictors", "eps", "minPts", "modelOpt", "seed",
                                                          "scaleEqualSD", "ready", "k-distplot"))
      jaspResults[["optimPlot"]] 		 $position <- 4
    }
  }
}





