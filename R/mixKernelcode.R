###   Jerôme MARIETTE - Nathalie VILLA-VIALANEIX code

#' @keywords internal
#' @importFrom quadprog solve.QP
#' @importFrom corrplot corrplot


umkl <- function(X, scale = TRUE, visualize=FALSE,
                 method = c("full", "statis", "sparse"), knn = 5,rho = 20) {

  #-- cosinus scaling --------------------#
  if (scale) {
    X.scaled <- lapply (X, function(x) {
      x.cosinus <- sweep(sweep(x, 2, sqrt(diag(x)), "/"),
                         1, sqrt(diag(x)), "/")
      t(t(x.cosinus - colSums(x.cosinus) / nrow(x.cosinus)) - rowSums(x.cosinus) /
          nrow(x.cosinus)) + sum(x.cosinus) / nrow(x.cosinus)^2
    })
  } else {
    X.scaled <- X
  }

  #-- UMKL approaches --------------------#
  beta <- 1 / length(X.scaled)
  if (method == 'statis') {

    similarities <- outer(1:length(X.scaled), 1:length(X.scaled),
                          FUN = Vectorize(function(i, j) {
                            sum(diag(X.scaled[[i]] %*% X.scaled[[j]])) / (norm(X.scaled[[i]], type="F") *
                                                                     norm(X.scaled[[j]], type="F"))
                          }))
    if(visualize) {
      rownames(similarities) <- colnames(similarities) <- names(X)
      corrplot(similarities, type = "full", tl.col = "black", tl.srt = 45,
               method = "pie")
      return(similarities)
    }
    weights <- eigen(similarities, symmetric = TRUE)$vectors[ ,1]
    weights <- weights / sum(weights)

  } else { # for both full-UMKL and sparse-UMKL methods

    # build the adjacency matrix considering the input kernels
    all.adjacency <- lapply(X.scaled, function(x) {
      adjacency.x <- apply(x, 1, function(row) {
        adjacency.row <- rank(row)
        adjacency.row[which(adjacency.row < length(adjacency.row)-knn)] <- 0
        adjacency.row[which(adjacency.row > 0)] <- 1
        adjacency.row
      })
      diag(adjacency.x) <- 0
      adjacency.x + t(adjacency.x)
    })
    graph.weights <- Reduce("+", all.adjacency)

    C.matrix <- outer(1:length(X.scaled), 1:length(X.scaled),
                      FUN = Vectorize(function(i, j) {
                        sum(graph.weights * apply(X.scaled[[i]], 1, function(rowi) {
                          apply(X.scaled[[i]], 1, function(rowj) {
                            sum(abs(rowi-rowj))
                          })
                        }) %*% apply(X.scaled[[j]], 1, function(rowi) {
                          apply(X.scaled[[j]], 1, function(rowj){
                            sum(abs(rowi-rowj))
                          })
                        }))
                      }))
    C.matrix <- matrix(unlist(C.matrix), nrow=length(X.scaled))

    if (method == 'sparse') { # if sparse-UMKL, use quadprog

      dvec <- rep(0, length(X.scaled))
      # cosinus normalization
      C.matrix.n <- sweep(sweep(C.matrix, 2, sqrt(diag(C.matrix)), "/"), 1,
                          sqrt(diag(C.matrix)), "/")
      # frobenius normalization
      # C.matrix.n <- C.matrix / norm(C.matrix, type="F")
      Dmat <- 2 * C.matrix.n
      Amat <- matrix(data = 1, nrow = nrow(Dmat), ncol = 1)
      Amat <- cbind(Amat, diag(1, nrow = nrow(Dmat), ncol = nrow(Dmat)))
      bvec <- rep(0, length(X.scaled) + 1)
      bvec[1] <- 1
      weights <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)$solution

    } else { # else, solve the problem using ADMM

      k <- 1
      Z <- rep(1/sqrt(length(X.scaled)), length(X.scaled))
      Y <- rep(0,length(X.scaled))
      threshold <- 10^(-2)
      # cosinus normalization
      C.matrix.n <- sweep(sweep(C.matrix, 2, sqrt(diag(C.matrix)), "/"), 1,
                          sqrt(diag(C.matrix)), "/")
      # frobenius normalization
      # C.matrix.n <- C.matrix / norm(C.matrix, type="F")
      Dmat <- C.matrix.n + diag(rho/2, ncol=length(X.scaled), nrow=length(X.scaled))
      Dmat <- 2 * Dmat
      Amat <- diag(1, nrow=dim(Dmat)[1], ncol=dim(Dmat)[1])
      bvec <- rep(0, length(X.scaled))
      repeat {
        # solve x using quadprog
        dvec <- rho*Z - Y
        X.ADMM <- solve.QP(Dmat, dvec, Amat, bvec, meq=0)$solution
        # project Z
        newZ <- X.ADMM / norm(matrix(X.ADMM), "2")
        if (k != 1 && norm(matrix(Z - newZ)) < threshold) break
        Z <- newZ
        # update Y
        Y <- Y + rho*(X.ADMM-Z)
        k <- k + 1
      }
      weights <- Z / norm(matrix(Z))

    }
  }
  return(weights)
}
