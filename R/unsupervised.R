### Unsupervised learning ###


#' Hierarchical clustering
#' @param data Input data: a matrix or data.frame.
#' @param kernel "linear" for linear kernel, "rbf" for RBF, "cRBF" for clrRBF, "qJac" for quantitative Jaccard.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param H Kernel hyperparameter.
#' @param method the agglomeration method to be used: "ward.D", "ward.D2" (Ward's (1963) clustering criterion),
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param plot TRUE to return the plot
#' @param labels A vector of labels for the leaves of the tree. Row names of
#' the input data are used by default. If set to FALSE, no labels at all are plotted.
#' @param title Plot title
#' @param cut an integer scalar or vector with the desired number of groups
#' @param height numeric scalar or vector with heights where the tree should be cut
#' @param colors border color(s) for the rectangles.
#' @return An object of class hclust with or without a plot of the cluster dendogram
#' @examples
#' hklust(soilDataRaw[-89,],kernel="cRBF", title = "Soil data cluster dendogram",cut=3,colors=2:4)
#' @importFrom stats as.dist hclust cutree rect.hclust
#' @export


hklust <- function(data, kernel="linear", H=NULL, method="ward.D2", plot=TRUE, labels=NULL, title=NULL, cut=NULL, height=NULL, colors="black") {
  kMatrix <- kernelSelect(data=data, kernel=kernel, h=H)
  distances <- as.dist(1-kMatrix)
  clustering <- hclust(distances,method = method)
  if(plot) {
   plot(clustering,main=title)
    if(!is.null(cut) || !is.null(height)) {
      clusters <- cutree(clustering, h = height, k = cut)
      rect.hclust(clustering , k = cut, h = height, border = colors)
    }
  }
  return(clustering)
}


#' kernel PCA
#' @param data Input data: a matrix or data.frame.
#' @param kernel "linear" for linear kernel, "rbf" for RBF, "cRBF" for clrRBF, "qJac" for quantitative Jaccard.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param H Kernel hyperparameter.
#' @param y Response vector; optional. Colors each dot (individual) according to its value in the response / target variable.
#' @param dim The two PC that have to be displayed. Defaults to the two first PC.
#' @param colors Dot fill color; optional.
#' @param na.col Depends on y and col. Dot fill color for the y missing values; defaults to grey.
#' @param title Plot title
#' @param legend Defaults to TRUE. A legend with the color corresponding to each y value.
#' @param labels If true, each dot will be labeled with its row number.
#' @return A  k-PCA plot generated with ggplot2.
#' @examples
#' kernPCA(soilDataRaw[-89,],kernel="cRBF", H=0.0001,y=soilMetaData$ph[-89], col=c("aquamarine2","orchid3"),
#' title = "Soil kernel PCA",legend = TRUE)
#' @importFrom catkern cplot
#' @export


kernPCA <- function(data, kernel, H=NULL, y, dim=c(1,2), colors, na.col="grey70", title, legend = TRUE, labels=FALSE) {
  kMatrix <- kernelSelect(data=data, kernel=kernel, h=H)
  return(cplot(matrix=kMatrix, y=y, dim=dim, col=colors, na.col=na.col, title=title, legend = legend, labels=labels))
}

