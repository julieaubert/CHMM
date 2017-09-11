# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Implementation of the Viterbi algorithm
#'
#'@param emisPr a matrix of emission probabilities for the considering series.
#'@param transPr a matrix of state transition probabilities.
#'@param initPr a vector specifying initial state probabilities.
#'@return \item{path}{the most likely path (state sequence).}
#'@export
#'@keywords internal

viterbi_algo <- function(emisPr, transPr, initPr){
  nb.states <- ncol(emisPr)
  nbT <- nrow(emisPr)
  logF <- matrix(0, nbT, nb.states)
  logF[1,] <- log(initPr) + log(emisPr[1,])

  for (tt in 2:nbT) for (k in 1:nb.states) logF[tt,k] <- log(emisPr[tt, k]) + max(logF[tt - 1,] + log(transPr[,k]))

  path <- rep(NA, nbT)
  path[nbT] <- which.max(logF[nbT,])
  for (tt in (nbT - 1):1) path[tt] <- which.max(logF[tt, ] + log(transPr[, path[tt + 1]]))
  return(path)
}

