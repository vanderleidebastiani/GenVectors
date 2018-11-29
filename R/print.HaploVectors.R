#' @rdname HaploVectors
#' @encoding UTF-8
#' @export
print.HaploVectors <- function(x, ...){
  cat("$call:\n")
  cat(deparse(x$call), "\n\n")
  cat("$haplotypes:\n")
  print(x$haplotypes, ...)
  cat("$haplotype.distances:\n")
  print(x$haplotype.distances, ...)
  cat("\n$haplotype.per.locality:\n")
  print(x$haplotype.per.locality, ...)
  cat("\n$haplotype.network:\n")
  print(x$haplotype.network, ...)
  cat("\n$vectors:\n")
  print(x$vectors, ...)
  cat("\n$values:\n")
  print(x$values, ...)
  cat("\n$correlations:\n")
  print(x$correlations, ...)
  cat("\n$P:\n")
  print(x$P, ...)
  cat("\n$scores.haplotypes:\n")
  print(x$scores.haplotypes, ...)
  if(!is.null(x$model)){
    cat("\n$model:\n")
    print(x$model, ...)
    cat("\n$obs.statistic:\n")
    print(x$obs.statistic, ...)
    cat("\n$p.site.shuffle:\n")
    print(x$p.site.shuffle, ...)
    cat("\n$p.taxa.shuffle:\n")
    print(x$p.taxa.shuffle, ...)  
  }
  invisible(x)
}
