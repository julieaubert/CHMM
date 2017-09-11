# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Perform variational inference of coupled Hidden Markov Models.
#'
#' @param X a data matrix of observations. Columns correspond to individuals.
#' @param nb.states a integer specifying the numbers of states.
#' @param S a matrix of similarity between individuals.
#' @param omega a value of omega.
#' @param meth.init a string specifying the initialization method ("mclust" or "kmeans"). The default method is "mclust".
#' @param var.equal a logical variable indicating whether to treat the variances as being equal.
#' @param itmax an integer specifying the maximal number of iterations for the EM algorithm.
#' @param threshold a value for the threshold used for the stopping criteria.
#' @return a list of 9 components
#' \describe{
#' \item{\code{postPr}}{a list containing for each series the posterior probabilities. }
#' \item{\code{initPr}}{ a numeric specifying the initial state probabilities.}
#' \item{\code{transPr}}{ a matrix of the state transition probabilities.}
#' \item{\code{esAvg}}{ a numeric of the estimated mean for each state.}
#' \item{\code{esVar}}{ a numeric of the estimated variance for each state.}
#' \item{\code{emisPr}}{a list containing for each series the emission probabilities.}
#' \item{\code{emisPrW}}{a list containing for each series the emission probabilities taking into account for the dependency structure.}
#' \item{\code{RSS}}{ a numeric corresponding to the Residuals Sum of Squares.}
#' \item{\code{iterstop}}{ an integer corresponding to the total number of iterations.}
#'}
#' @export
#' @references Wang, X., Lebarbier, E., Aubert, J. and Robin, S., Variational inference for coupled Hidden Markov Models applied to the joint detection of copy number variations.
#'

CHMM_VEM <- function(X, nb.states, S = NULL, omega = 0.7, meth.init = "mclust", var.equal = TRUE, itmax = 5e2,
                     threshold = 1e-7){
  if (is.null(S))
    stop("Argument S must be specified.")
  nbT <- nrow(X)
  nbI <- ncol(X)

  # initialisation   ------------------------------------------------------
 # if (is.null(init.esAvg)) {
    res.init <- init.VEM(X, nb.states, meth.init, var.equal, nbI, nbT)
    esAvg <- res.init$esAvg
    esVar <- res.init$esVar
    transPr <- res.init$transPr
    initPr <- res.init$initPr
    postPr.last <- res.init$postPr

  # Variational E-M algorithm    ------------------------------------------------------
  Old.param = c(esAvg, esVar, as.vector(transPr))
  for (iter in 1:itmax) {
    ## 1.VE-step    -----------------------------
    postPr.list <- list()
    emisPr.list <- list()
    emisPrW.list <- list()
    trsTmp <- matrix(0, nb.states, nb.states)
    emisPr <- matrix(0, nbT, nb.states)
    for (ind in 1:nbI) {
      w <- matrix(0, nbT, nb.states)
      for (ij in c(1:nbI)[-ind])
        w <- w + S[ind, ij] * (1 - postPr.last[[ij]])
      emisPr <- Emis.Gauss(X[,ind], esAvg, esVar)
      emisPr.list[[ind]] <- emisPr
      emisPrW <- omega^w * emisPr
      emisPrW.list[[ind]] <- omega^w * emisPr

      # Transform matrix to vectors     -----------------------------
      emisWVec <- as.vector(emisPrW)
      trsVec <- as.vector(transPr)
      initTmp <- as.vector(initPr * emisPrW[1,])
      initPr <- initTmp /sum(initTmp)
      # Forward-Backward recursion      -----------------------------
      resF <- ForwardR(emisWVec, initPr, trsVec)
      resB <- BackwardR(resF$Fpr, trsVec, nb.states)

      postPr.tmp <- matrix(resB$postPr, nbT, nb.states)
      postPr.tmp <- apply(postPr.tmp, 2, pmax, 1e-6)
      postPr <- postPr.tmp / rowSums(postPr.tmp)
      postPr.list[[ind]] <- postPr
      Fpr <- matrix(c(resF$Fpr), nbT, nb.states)
      Gpr <- matrix(c(resB$Gpr), nbT, nb.states)

      trsTmp <- trsTmp + transPr * t(Fpr[-nbT,]) %*% (postPr[-1,] / Gpr[-1,])
    }
    postPr.last <- postPr.list

    ## 2.M-step      ----------------------------
    # update transPr
    trsAvg <- trsTmp / nbI
    transPr <- trsAvg / rowSums(trsAvg)

    # update initPr       ----------------------------
    init.tmp <- rep(0, nb.states)
    for(i in 1:nbI) init.tmp <- init.tmp + postPr.list[[i]][1,]
    initPr <- init.tmp / sum(init.tmp)

    # update esAvg        ----------------------------
    esAvg <- NULL
    for (r in 1:nb.states) {
      nom <- 0
      den <- 0
      for (i in 1:nbI) {
        nom <- nom + sum(postPr.list[[i]][,r] * (X[,i]))
        den <- den + sum(postPr.list[[i]][,r])
        }
      esAvg <- c(esAvg, nom / den)
    }

    # update esVar         ----------------------------
    RSS <- 0
    esVar <- NULL
    if (var.equal == TRUE) {
      nom <- 0
      for(i in 1:nbI) nom <- nom + sum(X[,i]^2) - 2 * esAvg %*% t(postPr.list[[i]]) %*% X[,i] + sum(esAvg^2 %*% t(postPr.list[[i]]))
      esVar <- rep(nom/(nbI*nbT),nb.states)
      RSS <- nom
  #  } else if(var.model == "V") {
    } else {
      for (r in 1:nb.states) {
        nom <- 0
        den <- 0
        for (i in 1:nbI) {
          nom <- nom + sum(postPr.list[[i]][,r] * (X[,i]-esAvg[r])^2)
          den <- den + sum(postPr.list[[i]][,r])
          }
        esVar <- c(esVar, nom / den)
        RSS <- RSS + nom
      }
    }

    #  order the parameters          ----------------------------
    ordAvg <- order(esAvg)
    esAvg <- sort(esAvg)
    esVar <- esVar[ordAvg]
    transPr <- transPr[ordAvg,ordAvg]

    ## 3.Stop iteration          ----------------------------
    New.param <-  c(esAvg, esVar, as.vector(transPr))
    crit <- New.param - Old.param
    esVar <- pmax(esVar, 1e-3)
    Old.param  <- New.param

    if (iter > 1 && max(abs(crit)) <= threshold) break()
  }
  return(list(postPr = postPr.list, initPr = initPr,  transPr = transPr,  esAvg = esAvg, esVar = esVar,
              emisPr = emisPr.list, emisPrW = emisPrW.list, RSS = as.numeric(RSS), iterstop = iter))
}

