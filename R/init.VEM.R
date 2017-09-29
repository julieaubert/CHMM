# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
## __________________________________________________________
##
## Function :: init.VEM()
## __________________________________________________________
##
#' Initialization step of the \code{CHMM_VEM} function.
#'
#' @param X a matrix of observations. Columns correspond to series (individuals).
#' @param nb.states an integer specifying the numbers of states.
#' @param meth.init a string specifying the initialization method ("mclust" or "kmeans"). The default method is "mclust".
#' @param var.equal a logical variable indicating whether to treat the variances as being equal (TRUE, value by default) or not (FALSE).
#' @param nbI an integer specifying the number of series.
#' @param nbT an integer specifying the length of one series.
#'
#' @return A list containing the parameters of the model
#' \describe{
#' \item{\code{esAvg}}{ a numeric of the estimated mean for each state.}
#' \item{\code{esVar}}{ a numeric of the estimated variance for each state.}
#' \item{\code{transPr}}{ a matrix of the state transition probabilities}
#' \item{\code{postPr}}{ a list containing for each series the posterior probabilities.}
#' \item{\code{initPr}}{ a numeric specifying the initial state probabilities}.
#'}
#'
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom stats kmeans
#' @importFrom stats runif
#' @export




init.VEM <- function(X, nb.states, meth.init, var.equal, nbI, nbT){
  # mclust initialization ---------------------------------------------------------
  if(meth.init == "mclust") {
    if (var.equal == TRUE) {
      clust <- mclust::Mclust(data = as.vector(t(X)), G = nb.states, modelNames = "E")
      esVar <- rep(clust$parameters$variance$sigmasq, nb.states)
    } else {
      clust <- mclust::Mclust(data = as.vector(t(X)), G = nb.states, modelNames = "V")
      esVar <- clust$parameters$variance$sigmasq
    }
    esAvg <- clust$parameters$mean


    # kmeans initialization ---------------------------------------------------------
  } else if(meth.init == "kmeans") {
    clust <- kmeans(x = as.vector(t(X)), centers = nb.states)
    esAvg <- sort(as.vector(clust$centers))
    if (var.equal == TRUE) {
      esVar <- rep(clust$tot.withinss / (nbT * nbI - 1), nb.states)
    } else {
      esVar <- clust$withinss / (clust$size - 1)
    }
  }
  mat.tmp <- matrix(runif(nb.states^2), ncol = nb.states) + diag(rep(50, nb.states))
  transPr <- mat.tmp / rowSums(mat.tmp)

  # initial distribution -----------------------------------------
  eigenvalues <- round(eigen(t(transPr))$values, 3)
  pos  <- which(eigenvalues == 1.000)
  nuHMM <- eigen(t(transPr))$vectors[, pos]
  initPr <- pmax(as.numeric(nuHMM / sum(nuHMM)), 0)
  initPr <- initPr / sum(initPr)

  # postPr  ------------------------------------------------------
  postPr <- list()
  for(ind in 1:nbI){
    #tau.tmp <- data.frame(matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states))
    tau.tmp <- matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states)
    tau.tmp <- tau.tmp / rowSums(tau.tmp)
    postPr[[ind]] <- tau.tmp
  }
  return(list(esAvg = esAvg, esVar = esVar, transPr = transPr, initPr = initPr,
              postPr = postPr))

}
