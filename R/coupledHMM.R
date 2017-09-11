# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Perform inference of coupled hidden markov models.
#'
#'
#' @param X a matrix of observations. Columns correspond to series (individuals).
#' @param nb.states a integer specifying the numbers of states.
#' @param S a matrix of similarity between individuals.
#' @param omega.list a vector of omega values.
#' @param var.equal a logical variable indicating whether to treat the variances as being equal (var.equal = TRUE).
#' @param exact a logical variable indicating whether to use VEM (exact = FALSE) or EM (exact = TRUE) algorithm for the inference of the model.
#' @param meth.init a string specifying the initialization method ("mclust" or "kmeans") for the (V)-EM algorithm. The default method is "mclust".
#' @param viterbi a logical variable indicating whether to use Maximum A Posteriori method (FALSE) or Viterbi algorithm (TRUE, by default) for recovering the most likely path.
#' @param itmax an integer specifying the maximal number of iterations for the CHMM_(V)EM algorithm.
#' @param threshold a value for the threshold used for the stopping criteria for the CHMM_(V)EM algorithm.
#' @return A list of 4 objets.
#' \describe{
#' \item{\code{omega}}{ an integer corresponding to the selected value among the omega.list.}
#' \item{\code{model}}{a list corresponding to the output of the \code{CHMM-EM} or \code{CHMM-VEM} function for the selected model.}
#' \item{\code{status}}{a matrix with status associated to each series in column and each position in row.}
#' \item{\code{RSS.omega}}{ a dataframe with omega values and the associated Residuals Sum of Squares.}
#'}
#' @export
#' @examples
#' data(toyexample)
#' # Variational inference of a coupled hidden Markov Chains
#' resCHMM <- coupledHMM(X = toydata, nb.states = 3, S = cor(toystatus),
#'                       omega.list = c(0.3, 0.5, 0.7, 0.9))
#' # Breakpoints positions and status of segments
#' info <- clusterseg(resCHMM$status)
#' # head(info)
#' @seealso \code{\link{CHMM_VEM}}, \code{\link{CHMM_EM}}
#' @references Wang, X., Lebarbier, E., Aubert, J. and Robin, S., Variational inference for coupled Hidden Markov Models applied to the joint detection of copy number variations.

coupledHMM <- function(X, nb.states = 3, S = NULL, omega.list = c(0.3, 0.7, 0.9),
         var.equal = TRUE, exact = FALSE, meth.init = "mclust", viterbi = TRUE, itmax = 5e2,
         threshold = 1e-7){
  if (is.null(S))
      stop("Argument S must be specified.")
  nbI <- ncol(X)
  # inference ---------------------------
  RSS <- Inf
  if(exact == FALSE){
    resApprox <- apply(as.matrix(omega.list),c(1),FUN=function(y) CHMM_VEM(X, nb.states, S , y, meth.init, var.equal, itmax, threshold))
    RSS <- unlist(lapply(resApprox,FUN = function(x) x$'RSS'))
    mod.sel <- resApprox[[which.min(RSS)]]
    omega.sel <- omega.list[which.min(RSS)]

    if (viterbi == TRUE){
      status <- apply(as.matrix(1:nbI),c(1),FUN=function(y) viterbi_algo(mod.sel$emisPrW[[y]], mod.sel$transPr, mod.sel$initPr))
    }else{
      status <-  apply(as.matrix(1:nbI),c(1),FUN=function(y) apply(mod.sel$postPr[[y]], 1, which.max))
    }
    colnames(status) <- colnames(X)
  }else{

    resExact <- apply(as.matrix(omega.list),c(1),FUN=function(y)
      CHMM_EM(X, nb.states, S , y, meth.init, var.equal, itmax, threshold))
    RSS <- unlist(lapply(resExact,FUN = function(x) x$'RSS'))
    mod.sel <- resExact[[which.min(RSS)]]
    omega.sel <- omega.list[which.min(RSS)]

    if (viterbi == TRUE){
      stsGb <- viterbi_algo(mod.sel$'emisGb', mod.sel$'transGb', mod.sel$'initGb')
      status <- mod.sel$ID.K[stsGb, ]
    }else{
      status <-  apply(as.matrix(1:nbI),c(1),FUN=function(y) apply(mod.sel$postPr[[y]], 1, which.max))
    }

    colnames(status) <- colnames(X)
  }#end if
  return(list(omega = omega.sel, model = mod.sel, status = status, RSS.omega = data.frame(RSS = RSS, omega = omega.list)))
}
