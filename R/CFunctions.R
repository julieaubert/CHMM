# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
# __________________________________________________________
#
# ForwardR
# __________________________________________________________
#' Forward step
#'
#' @param emisVec a vector of emission probabilities.
#' @param initPr a vector specifying initial state probabilities.
#' @param trsVec a vector of state transition probabilities.
#' @useDynLib CHMM Forward
#' @export
#' @keywords internal
ForwardR <- function(emisVec,initPr,trsVec) {
  nb.states <- length(initPr)
    nbI <- length(emisVec) / nb.states
    .C("Forward", as.double(emisVec), as.double(initPr), as.double(trsVec),
      as.integer(nbI), as.integer(nb.states), Fpr = double(nbI * nb.states),
       Lambda = double(nbI), PACKAGE = "CHMM")
}
## __________________________________________________________
##
## BackwardR
## __________________________________________________________
##
#' Backward step
#'
#' @param Fpr Fpr.
#' @param trsVec a vector of state transition probabilities.
#' @param nb.states an integer specifying the numbers of states.
#' @useDynLib CHMM Backward
#' @export
#' @keywords internal
BackwardR <- function(Fpr, trsVec, nb.states) {
    nbI <- length(Fpr) / nb.states
    .C("Backward", as.double(Fpr), as.double(trsVec), as.integer(nbI),
       as.integer(nb.states), postPr = double(nbI * nb.states), Gpr = double(nbI * nb.states), PACKAGE = "CHMM")
    }

