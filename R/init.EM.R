# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
## __________________________________________________________
##
## Function :: init.EM()
## __________________________________________________________
##
#' Initialization step of the \code{CHMM_EM} function.
#'
#' @param X a matrix of observations. Columns correspond to series (individuals).
#' @param nb.states an integer specifying the numbers of states.
#' @param meth.init a string specifying the initialization method ("mclust" or "kmeans"). The default method is "mclust".
#' @param var.equal a logical variable indicating whether to treat the variances as being equal (TRUE, value by default) or not (FALSE).
#' @param nbI an integer specifying the number of series.
#' @param nbT an integer specifying the length of one series.
#' @return A list of 6 objects.
#' \describe{
#' \item{\code{esAvgGb}}{ a matrix of \code{nbK}(nb.states^nbI) rows and \code{nbI} columns of estimated mean.}
#' \item{\code{esVarGb}}{ a matrix of \code{nbK}(nb.states^nbI) rows and \code{nbI} columns of estimated variance.}
#' \item{\code{esAvg}}{ a numeric of the estimated mean for each state.}
#' \item{\code{esVar}}{ a numeric of the estimated variance for each state.}
#' \item{\code{transGb}}{ a matrix of the state transition probabilities.}
#' \item{\code{initGb}}{ a numeric specifying the initial state probabilities.}
#'}
#'@details By default, an initialization with the \code{meth.init="mclust"} is performed with homogeneous variances.
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom stats kmeans
#' @importFrom stats runif
#' @export
#' @seealso \code{\link{CHMM_EM}}
#'
init.EM <- function(X, nb.states, meth.init, var.equal, nbI, nbT){
  nbK <- nb.states^nbI
  ID.K <- ExpGrid(nb.states, nbI)
  # mclust initialization ---------------------------------------------------------
  if(meth.init == "mclust") {
    if (var.equal == TRUE) {
      clust <- mclust::Mclust(data = as.vector(t(X)), G = nb.states, modelNames = "E")
      esVar <- clust$parameters$variance$sigmasq
      esVarGb <- matrix(esVar, nbK, nbI)
    } else {
      clust <- mclust::Mclust(data = as.vector(t(X)), G = nb.states, modelNames = "V")
      esVar <- clust$parameters$variance$sigmasq
      esVarGb <- matrix(esVar[ID.K], nbK, nbI)
    }
    esAvg <- clust$parameters$mean
# kmeans initialization ---------------------------------------------------------
   } else if(meth.init == "kmeans") {
    clust <- kmeans(x = as.vector(t(X)), centers = nb.states)
    esAvg <- sort(as.vector(clust$centers))
    if (var.equal == TRUE) {
      esVar <- clust$tot.withinss / (nbT * nbI - 1)
      esVarGb <- matrix(esVar, nbK, nbI)
    } else {
      esVar <- clust$withinss / (clust$size - 1)
      esVarGb <- matrix(esVar[ID.K], nbK, nbI)
    }
  }
  esAvgGb <- matrix(esAvg[ID.K], nbK, nbI)

  # transGb -----------------------------------------------------------
  mat.tmp  <- matrix(runif(nbK^2), nbK, nbK) + diag(rep(50, nbK))
  transGb <- mat.tmp / rowSums(mat.tmp)

  ## initGb -----------------------------------------------------------
  val.propre  <- round(eigen(t(transGb))$values, 3)
  pos  <- which(val.propre == 1.000)
  nuHMM  <- eigen(t(transGb))$vectors[,pos]
  nuHMM  <- nuHMM / sum(nuHMM)
  initGb  <- pmax(as.numeric(nuHMM), 0)
  initGb <- initGb / sum(initGb)

   return(list(esAvgGb = esAvgGb, esVarGb = esVarGb, esAvg = esAvg, esVar = esVar, transGb = transGb, initGb = initGb))

}


