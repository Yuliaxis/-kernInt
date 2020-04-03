### Unsupervised learning ###


#' Hierarchical clustering
#' @param data Input data
#' @param kernel "linear" for linear kernel, "rbf" for RBF, "cRBF" for clrRBF, "qJac" for quantitative Jaccard.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param comb If data is a list or array: how to combine them ("mean","statis","full","sparse" or a vector of coefficients)
#' @param H Kernel gamma hyperparameter if needed (only RBF-like kernels)
#' @param domain Only used in "frbf" or "flin".
#' @param method the agglomeration method to be used: "ward.D", "ward.D2" (Ward's (1963) clustering criterion),
#' "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param plot TRUE to return the plot, FALSE to return the dendogram object
#' @param labels A vector of labels for the leaves of the tree. Row names of
#' the input data are used by default. If set to FALSE, no labels at all are plotted.
#' @param title Plot title
#' @param cut an integer scalar or vector with the desired number of groups
#' @param height numeric scalar or vector with heights where the tree should be cut
#' @param colors border color(s) for the rectangles.
#' @return An object of class hclust with or without a plot of the cluster dendogram
#' @examples
#' hklust(soil$abund,kernel="jac", title = "Soil data cluster dendogram",cut=3,colors=2:4)
#' @importFrom stats as.dist hclust cutree rect.hclust
#' @export


hklust <- function(data, comb="mean", kernel, H=NULL, domain=NULL, method="ward.D2", plot=TRUE, labels=NULL, title=NULL, cut=NULL, height=NULL, colors="black") {
  if(class(data) == "list" | class(data) == "array" ) {
    kMatrix <-fuseData(DATA=data,kernels=kernel,coeff=comb,h=H)
  }  else  {
   kMatrix <- kernelSelect(data=data, kernel=kernel, h=H)
  }
  distances <- stats::as.dist(1-kMatrix)
  clustering <- stats::hclust(distances,method = method)
  if(plot) {
   plot(clustering,main=title)
    if(!is.null(cut) || !is.null(height)) {
      clusters <- stats::cutree(clustering, h = height, k = cut)
      stats::rect.hclust(clustering , k = cut, h = height, border = colors)
    }
  }
  return(clustering)
}


#' kernel PCA
#' @param data Input data: a matrix or data.frame.
#' @param kernel "linear" for linear kernel, "rbf" for RBF, "cRBF" for clrRBF, "qJac" for quantitative Jaccard.
#' "matrix" if a pre-calculated kernel matrix is given as input.
#' @param plot TRUE to return the plot, FALSE to return the projection object
#' @param H Kernel gamma hyperparameter if needed (only RBF-like kernels)
#' @param domain Only used in "frbf" or "flin".
#' @param y Response vector; optional. Colors each dot (individual) according to its value in the response / target variable.
#' @param dim The two PC that have to be displayed. Defaults to the two first PC.
#' @param colors Dot fill color; optional.
#' @param na.col Depends on y and colors. Dot fill color for the y missing values; defaults to grey.
#' @param title Plot title
#' @param legend Defaults to TRUE. A legend with the color corresponding to each y value.
#' @param labels If true, each dot will be labeled with its row number.
#' @return A  k-PCA plot generated with ggplot2.
#' @examples
#' kernPCA(soil$abund,kernel="clin", y=soil$metadata$phd, colors=c("aquamarine2","orchid3"),
#' title = "Soil kernel PCA",legend = TRUE)
#' @import ggplot2
#' @importFrom  graphics plot
#' @export


kernPCA <- function(data, kernel, plot=TRUE,H=NULL,domain=NULL, y, dim=c(1,2), colors, na.col="grey70", title, legend = TRUE, labels=FALSE) {
  matrix <- kernelSelect(data=data, kernel=kernel, h=H)
  matrix <- kernlab::as.kernelMatrix(matrix)
  i <- dim[1]; j <- dim[2]
  subjects.kpca <- kernlab::kpca(matrix,  kernel = matrix)
  xpercent <- kernlab::eig(subjects.kpca)[i]/sum(kernlab::eig(subjects.kpca))*100
  ypercent <- kernlab::eig(subjects.kpca)[j]/sum(kernlab::eig(subjects.kpca))*100
  df <- as.data.frame(kernlab::rotated(subjects.kpca))
  if(!plot) return(df)
  q <- ggplot2::ggplot(df,ggplot2::aes(df[,i], df[,j])) + ggplot2::theme_bw() +
    ggplot2::xlab(paste("PC", i," - ", format(xpercent,digits=3), "%",sep="")) + ggplot2::ylab(paste("PC", j," - ", format(ypercent,digits=3), "%",sep=""))

  if(hasArg(y)) {
    q <- q + ggplot2::geom_point(ggplot2::aes(colour = y))
    if(hasArg(colors)) q <- q + ggplot2::scale_colour_gradientn(colours=colors, na.value = na.col)
    if(!legend) q <- q + ggplot2::theme(legend.position = "none")
  } else if(hasArg(colors)) {
    q <-  q + ggplot2::geom_point(colour = colors)
  } else {
    q <-  q + ggplot2::geom_point()
  }
  if(hasArg(title)) q <-  q + ggplot2::ggtitle(title)
  if(labels) q <- q + ggplot2::geom_text(ggplot2::aes(label=row.names(matrix),hjust=0, vjust=0))
  return(q)
}

