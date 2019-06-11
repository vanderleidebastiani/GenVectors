#' @rdname GenVectors
#' @encoding UTF-8
#' @export
print.GenVectors <- function(x, ...){
  cat("$call:\n")
  cat(deparse(x$call), "\n\n")
  if(!is.null(x$haplotypes)){
    cat("$haplotypes:\n")
    print(x$haplotypes, ...)
    cat("$haplotype.distances:\n")
    print(x$haplotype.distances, ...)
    cat("\n$individual.per.haplotype:\n")
    print(x$individual.per.haplotype, ...)
    cat("\n$haplotype.per.locality:\n")
    print(x$haplotype.per.locality, ...)
  }
  if(!is.null(x$SNP.distances)){
    cat("$SNP.distances:\n")
    print(x$SNP.distances, ...)
    cat("$pop:\n")
    print(x$pop, ...)
  }
  if(!is.null(x$genetic.distances)){
    cat("$genetic.distances:\n")
    print(x$genetic.distances, ...)
    cat("$pop:\n")
    print(x$pop, ...)
  }
  cat("\n$vectors:\n")
  print(x$vectors, ...)
  cat("\n$values:\n")
  print(x$values, ...)
  cat("\n$correlations:\n")
  print(x$correlations, ...)
  cat("\n$P:\n")
  print(x$P, ...)
  cat("\n$scores:\n")
  print(x$scores, ...)
  if(!is.null(x$model)){
    cat("\n$model:\n")
    print(x$model, ...)
    cat("\n$obs.statistic:\n")
    print(x$obs.statistic, ...)
    cat("\n$p.turnover:\n")
    print(x$p.turnover, ...)
    cat("\n$p.divergence:\n")
    print(x$p.divergence, ...)  
  }
  invisible(x)
}