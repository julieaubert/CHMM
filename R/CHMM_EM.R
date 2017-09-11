# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Perform exact inference of coupled hidden markov models.
#'
#' @param X a data matrix of observations. Columns correspond to individuals.
#' @param nb.states a integer specifying the numbers of states.
#' @param S a matrix of similarity between individuals.
#' @param omega a value of omega.
#' @param meth.init a string specifying the initialization method ("mclust" or "kmeans"). The default method is "mclust".
#' @param var.equal a logical variable indicating whether to treat the variances as being equal.
#' @param itmax an integer specifying the maximal number of iterations for the EM algorithm.
#' @param threshold a value for the threshold used for the stopping criteria.
#' @return a list of 10 components
#' \describe{
#' \item{\code{postPr}}{a list containing for each series the posterior probabilities.}
#' \item{\code{initGb}}{a numeric specifying the initial state probabilities.}
#' \item{\code{transGb}}{a matrix of the state transition probabilities.}
#' \item{\code{emisGb}}{a list containing for each series the emission probabilities.}
#' \item{\code{esAvg}}{ a numeric of the estimated mean for each state.}
#' \item{\code{esVar}}{ a numeric of the estimated variance for each state.}
#' \item{\code{ID.K}}{ a matrix containing all combination of possible state for nbI series.}
#' \item{\code{loglik}}{ a numeric with the value of the loglikelihood.}
#' \item{\code{RSS}}{ a numeric corresponding to the Residuals Sum of Squares.}
#' \item{\code{iterstop}}{ an integer corresponding to the total number of iterations.}
#'}
#' @export
#' @references Wang, X., Lebarbier, E., Aubert, J. and Robin, S., Variational inference for coupled Hidden Markov Models applied to the joint detection of copy number variations.
#'


## __________________________________________________________
##
## Function :: CHMM.EM()
## __________________________________________________________
##

CHMM_EM <- function(X, nb.states, S, omega, meth.init = "mclust", var.equal = TRUE,
                    itmax = 5e2, threshold = 1e-7){
    if (is.null(S))
    stop("Argument S must be specified.")
    nbT <- nrow(X)
    nbI <- ncol(X)
    nbK <- nb.states^nbI
    ID.K <- ExpGrid(nb.states, nbI)
    wPr <- Wl(ID.K, S, omega)

    # Initialization EM ---------------------------------------------
    res.init <- init.EM(X, nb.states, meth.init, var.equal, nbI, nbT)
    esAvgGb <- res.init$esAvgGb
    esVarGb <- res.init$esVarGb
    esVar <- res.init$esVar
    esAvg <- res.init$esAvg
    transGb <- res.init$transGb
    initGb <- res.init$initGb

    # E-M algorithm  ------------------------------------------------
    # Attention: besoin de esVar et esAvg
    Old.param <- c(esVar, esAvg)
    for (iter in 1:itmax) {
      emisGb <- EmisGb.Gauss(X, esAvgGb, esVarGb)
      # E-Step -------------------------------------------------
      ## Transform matrix to vectors  ---------------------------
      trsGbVec <- as.vector(transGb)
      # Forward-Backward recursion   ---------------------------
      resF <- ForwardR(as.vector(emisGb), initGb, trsGbVec)
      resB <- BackwardR(resF$Fpr, trsGbVec, nbK)

      postGb.tmp <- apply(matrix(resB$postPr, nbT, nbK), 2, pmax, 1e-6)
      postGb <- postGb.tmp / rowSums(postGb.tmp)
      Fpr <- matrix(c(resF$Fpr), nbT, nbK)
      Gpr <- matrix(c(resB$Gpr), nbT, nbK)
      loglik <- -sum(log(resF$Lambda))

      ##  M-Step ---------------------------
      ### update global mean: esAvgGb  ---------------------------
      for (q in 1:nb.states) {
        tauq <- NULL
        for (ind in 1:nbI) {
          idx <- which(ID.K[,ind] == q)
          tauq <- cbind(tauq, rowSums(postGb[,idx]))
          }
        esAvg[q] <- sum(tauq * X) / sum(tauq)
      }
      esAvg <- sort(esAvg)
      esAvgGb <- matrix(esAvg[ID.K], nbK, nbI)

      ### update global variance: esVarGb   ---------------------------
      if (var.equal == TRUE) {
        var.tmp <- 0
        for(q in 1:nb.states){
          tauiq <- NULL
          for (ind in 1:nbI) {
            idx <- which(ID.K[, ind] == q)
            tauiq <- cbind(tauiq, rowSums(postGb[, idx]))
            }
          var.tmp <- var.tmp + sum(tauiq * (X - esAvg[q]) ** 2)
        }
        esVar <- var.tmp / nbI / nbT
        esVarGb <- matrix(esVar, nbK, nbI)
      } else {
        for(q in 1:nb.states){
          tauiq <- NULL
          for (ind in 1:nbI) {
            idx <- which(ID.K[, ind] == q)
            tauiq <- cbind(tauiq, rowSums(postGb[,idx]))
            }
          esVar[q] <- sum(tauiq * ( X - esAvg[q]) ** 2) / sum(tauiq)
        }
        esVarGb <- matrix(esVar[ID.K], nbK, nbI)
      }

      ### update transition    ---------------------------
      Nkl <- transGb * (t(Fpr[-nbT,]) %*% (postGb[-1,] / Gpr[-1,]))
      Pkl <- t(Nkl) / wPr
      Pkl <- t(Pkl)

      nqr <- matrix(0, nb.states, nb.states)
      for (iq in 1:nb.states) {
        for (ir in 1:nb.states) {
          Ntld <- matrix(0, nbK, nbK)
          for (ik in 1:nbK) {
            for (il in 1:nbK) {
              tmp <- 0
              for (ii in 1:nbI) {
                qik <- ID.K[ik,ii]
                qil <- ID.K[il,ii]
                if (qik == iq && qil == ir) tmp = tmp + 1
                }
              Ntld[ik,il] <- Pkl[ik,il] * tmp
            }
          }
          nqr[iq,ir] <- sum(Ntld)
        }
      }
      transPr <- nqr / rowSums(nqr)

      trans.tmp <- matrix(0,nbK,nbK)
      for (k in 1:nbK) {
        qvec <- ID.K[k,]
        for (l in 1:nbK) {
          rvec <- ID.K[l,]
          trans.tmp[k,l] <- prod(diag(transPr[qvec,rvec])) * wPr[l]
        }
        }
      transGb <- trans.tmp / rowSums(trans.tmp)

      ### update initial distribution     ---------------------------
      val.propre <- round(eigen(t(transGb))$values, 3)
      pos <- which(val.propre == 1.000)
      muHMM <- eigen(t(transGb))$vectors[,pos]
      muHMM <- muHMM / sum(muHMM)
      init.tmp <- as.numeric(muHMM)
      init.tmp[init.tmp < 0] <- 0
      initGb  <- init.tmp / sum(init.tmp)

      ## Stop criteria      ---------------------------
      New.param <- c(esVar, esAvg)
      crit <- New.param - Old.param
      # esVarGb      = pmax(esVarGb,1e-3)
      Old.param <- New.param

      if (iter > 1 && max(abs(crit)) <= threshold) break()
    }#end for
    postPr.list = list()
    for (i in 1:nbI) {
      postTmp <- NULL
      for (q in 1:nb.states) {
        idq <- which(ID.K[,i] == q)
        postTmp <- cbind(postTmp, rowSums(postGb[,idq]))
        }
      postPr.list[[i]] <- postTmp
    }
    RSS <- 0
    if (var.equal == TRUE) {
      for (i in 1:nbI)
        RSS <- RSS + sum(X[,i]^2) - 2 * esAvg %*% t(postPr.list[[i]]) %*% X[,i]
      + sum(esAvg^2 %*% t(postPr.list[[i]]))
    } else {
      for (r in 1:nb.states)
        for(i in 1:nbI)
          RSS <- RSS + sum(postPr.list[[i]][,r] * (X[,i] - esAvg[r])^2)
    }
    return(list(postPr = postPr.list, initGb = initGb, transGb = transGb,
                emisGb = emisGb, esAvg = esAvg, esVar = esVar,
                ID.K = ID.K, loglik = loglik, RSS = as.numeric(RSS), iterstop = iter))
}












