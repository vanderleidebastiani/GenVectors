#' @rdname GenVectors
#' @encoding UTF-8
#' @export
SNPVectors <- function(x, pop, checkdata = TRUE, method = "euclidean", squareroot.dis = TRUE, choices = c(1, 2), analysis = "none", analysis.method = "none", envir, choices.analysis, runs = 999, ...){
  res <- SNPDist(x)
  res$call <- match.call()
  if (checkdata) {
    if (is.null(colnames(pop))) {
      stop("\n Erro in column names of pop data\n")
    }
    match.names <- match(colnames(res$SNP.distances), colnames(pop))
    if (sum(is.na(match.names)) > 0) {
      print("There are individuals from x data that are not on pop:", quote = FALSE)
      print(setdiff(colnames(res$SNP.distances), colnames(pop)))
      stop("\n Individuals not found on pop matrix\n")
    }
    pop <- as.matrix(pop[, match.names, drop = FALSE])
  }
  res$pop <- pop
  res.eigen <- PCPS::pcps(res$pop, phylodist = res$SNP.distances, method = method, squareroot = squareroot.dis)
  res.eigen$call <- NULL
  rownames(res.eigen$values) <- sub("pcps", "SNPvector", rownames(res.eigen$values))
  colnames(res.eigen$vectors) <- sub("pcps", "SNPvector", colnames(res.eigen$vectors))
  colnames(res.eigen$correlations) <- sub("pcps", "SNPvector", colnames(res.eigen$correlations))
  res <- c(res, res.eigen)
  res$scores <- summary(res.eigen, choices = choices)$scores$scores.species
  # Analysis <- c("none", "SNPmatrix", "SNPvector")
  Analysis <- c("none", "matrix", "vector")
  analysis <- pmatch(analysis, Analysis)
  if (length(analysis) != 1 | (is.na(analysis[1]))) {
    stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
  }
  FUN <- PCPS::select.pcpsmethod(analysis.method)
  if(analysis!=1 & !is.null(FUN)){
    if(analysis == 2){
      test <- PCPS::matrix.p.sig(pop, phylodist = res$SNP.distances, method.p = method, sqrt.p = squareroot.dis, FUN = FUN, envir = envir, runs = runs, newname = "SNPvector", ...)
    }
    if(analysis == 3){
      test <- PCPS::pcps.sig(pop, phylodist = res$SNP.distances, method = method, squareroot = squareroot.dis, choices = choices.analysis, FUN = FUN, envir = envir, runs = runs, newname = "SNPvector", ...)
    }
    test$call <- NULL
    test$PCPS.obs <- NULL
    res <- c(res, test)
  }
  class(res) <- "GenVectors"
  return(res)
}