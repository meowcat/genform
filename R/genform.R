#' @useDynLib genform
#' @importFrom Rcpp sourceCpp
NULL

#' @export
genform <- function(ms, msms = NULL, m = NULL, settings = NULL) {
  GenFormMatchIsotopeMsMs_R(ms, msms, m, settings)
}

