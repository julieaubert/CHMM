# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Wl Internal function calculating Wl.
#'
#' @param nb.states an integer specifying the numbers of states.
#' @param nbI an integer specifying the number of series.
#' @return a matrix containing all combination of possible state for nbI series.
#' @export
#' @keywords internal
ExpGrid <- function(nb.states, nbI){
  KI.list <- list()
  for (i in 1:nbI) {
    KI.list[[i]] <- 1:nb.states
  }
  KI.grid <- as.matrix(expand.grid(KI.list))
  return(KI.grid[, nbI:1])
}
#' Wl Internal function calculating Wl.
#'
#' @param ID.K as.matrix(expand.grid(list(c(1:nbI),c(1:nb.states))).
#' @param S a matrix of similarities between individuals.
#' @param omega .
#' @return Wl.
#' @importFrom stats dist
#' @export
#' @keywords internal

Wl <- function(ID.K, S, omega){
  nbK <- nrow(ID.K)
  Sl <- sapply(1:nbK, function(i) sum(S[lower.tri(S)][which(dist(ID.K[i,])>0)]))
  Wl <- omega^Sl
  return(Wl)
}

#' Summarize the results of the coupled HMM.
#'
#' @param x a matrix of status. Columns corresponds to series (individuals).
#' @return a data.frame with 4 columns
#'  \describe{
#' \item{\code{sample}}{name of the sample (series). }
#' \item{\code{posbegin}}{beginning position.}
#' \item{\code{posend}}{ending position.}
#' \item{\code{status}}{status.}
#'}
#' @export
#'
clusterseg <- function(x){
  output <- NULL
  for (i in 1:ncol(x)){
    succStatus <- diff(x)
    posend <- which(succStatus[,i] != 0)
    posbegin <- c(1, posend + 1)
    tmp <- data.frame(sample = colnames(x)[i], posbegin, posend = c(which(succStatus[,i]!=0), nrow(x)), status = x[posbegin,i])
    output <- rbind(output,tmp)
  }
  return(output)
}


