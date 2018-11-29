#' @title Internal functions
#' 
#' @description Internal functions with small modifications of \code{\link{matrix.p.sig}} and \code{\link{pcps.sig}}.
#' 
#' @encoding UTF-8
#' @import PCPS
#' @importFrom parallel makeCluster stopCluster parLapply
#' @include matrix.p.sig.mod.R
#' @include pcps.sig.mod.R
#' @param x An object returned by \code{\link{matrix.p.null}}.
#' @param replacement A replacement name to matched in object returned by \code{\link{matrix.p.null}}.
#' @param newname New name to be replaced in object returned by \code{\link{matrix.p.null}}.
#' @param comm Community data, with species as columns and sampling units as rows.
#' @param phylodist Matrix containing phylogenetic distances between species.
#' @param method Dissimilarity index.
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index.
#' @param FUN An object of class function to perform the analysis. 
#' @param choices Numeric vector to choose the axis used in analysis.
#' @param runs Number of permutations for assessing significance.
#' @param parallel Number of parallel processes or a predefined socket cluster done with parallel package.
#' @param ... Other arguments passed to FUN function.
#' @export
mutate.names.matrix.p.null <- function(x, replacement, newname){
  f.mut <- function(x, replacement, newname){
    colnames(x) <- sub(replacement, newname, colnames(x))
    return(x)
  }
  if(!is.null(x$pcps.obs)){
    x$pcps.obs <- f.mut(x$pcps.obs, replacement, newname)
  }
  if(!is.null(x$pcps.null.site)){
    x$pcps.null.site <- sapply(x$pcps.null.site, f.mut, replacement = replacement, newname = newname, simplify = FALSE)
  }
  if(!is.null(x$pcps.null.taxa)){
    x$pcps.null.taxa <- sapply(x$pcps.null.taxa, f.mut, replacement = replacement, newname = newname, simplify = FALSE)
  }
  if(!is.null(x$pcps.null.taxa.adj)){
    x$pcps.null.taxa.adj <- sapply(x$pcps.null.taxa.adj, f.mut, replacement = replacement, newname = newname, simplify = FALSE)
  }
  return(x)
}