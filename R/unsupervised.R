### Unsupervised learning ###


#' Hierarchical clustering
#' @param data Input data
#' @param kernel  "lin" or rbf" to standard Linear and RBF kernels. "clin" for compositional linear and "crbf" for Aitchison-RBF
#' kernels. "jac" for quantitative Jaccard / Ruzicka kernel. "jsk" for Jensen-Shannon Kernel. "flin" and "frbf" for functional linear
#' and functional RBF kernels. "matrix" if a pre-computed kernel matrix is given as input.
#' With an array or a list of length *m*: Vector of *m* kernels to apply to each dataset.
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
#' ## Simple case
#' hklust(soil$abund,kernel="clin", title = "Soil data cluster dendogram",labels=TRUE)
#' ## Spatial data fusion case
#' hklust(data=smoker$abund,kernel=rep("clin",4),title="Nose samples",cut=2,colors=c("black","red"))
#' @importFrom stats as.dist hclust cutree rect.hclust
#' @importFrom methods is
#' @export


hklust <- function(data, comb="mean", kernel, H=NULL, domain=NULL, method="ward.D2", plot=TRUE, labels=NULL, title=NULL, cut=NULL, height=NULL, colors="black") {
  if(is(data,"data.frame") | is(data, "matrix") ) {
    kMatrix <- kernelSelect(data=data, kernel=kernel, h=H)
    }  else  {
    kMatrix <-fuseData(DATA=data,kernels=kernel,coeff=comb,h=H)
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
#' @param kernel lin" or rbf" to standard Linear and RBF kernels. "clin" for compositional linear and "crbf" for Aitchison-RBF
#' kernels. "jac" for quantitative Jaccard / Ruzicka kernel. "jsk" for Jensen-Shannon Kernel. "flin" and "frbf" for functional linear
#' and functional RBF kernels. "matrix" if a pre-computed kernel matrix is given as input.
#' With an array or a list of length *m*: Vector of *m* kernels to apply to each dataset.
#' @param comb If data is a list or array: how to combine them ("mean","statis","full","sparse" or a vector of coefficients)
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
#' kernPCA(soil$abund,kernel="clin", y=soil$metadata$phd, colors=c("orange","orchid3"),
#' title = "Soil kernel PCA",legend = TRUE)
#' ## Heterogeneous data fusion case
#' Airway <- list()
#' Airway$nosel <- CSSnorm(smoker$abund$oroL)
#' Airway$throatl <- CSSnorm(smoker$abund$oroR)
#' smoking <- smoker$metadata$smoker[seq(from=1,to=62*4,by=4)]
#' kernPCA(data=Airway,kernel=rep("jac",2),title="Airway samples",y=smoking)#'
#' @import ggplot2
#' @importFrom  graphics plot
#' @importFrom methods is

#' @export


kernPCA <- function(data, comb="mean",kernel, plot=TRUE,H=NULL,domain=NULL, y, dim=c(1,2), colors, na.col="grey70", title, legend = TRUE, labels=FALSE) {
  if(is(data,"data.frame") | is(data, "matrix") ) {
    kMatrix <- kernelSelect(data=data, kernel=kernel, h=H)
  }  else  {
    kMatrix <-fuseData(DATA=data,kernels=kernel,coeff=comb,h=H)
  }
  if(sum(grepl("rbf",kernel)>0) && is.null(H)) stop("H is mandatory")
  matrix <- kernlab::as.kernelMatrix(kMatrix)
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

