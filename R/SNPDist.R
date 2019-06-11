#' @title SNPDist 
#' 
#' @description Function to computes distances between alleles. This function use the function \code{\link{propShared}} of package adegenet.
#' 
#' @encoding UTF-8
#' @importFrom adegenet propShared
#' @param x A set of individual genotypes (an object of class \code{\link{genind}}) as used by the function \code{\link{propShared}}.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{SNP.distances}{A matrix of distances between alleles.} 
#' @seealso \code{\link{SNPVectors}}
#' @export
SNPDist <- function(x){
  res <- list(call = match.call())
  res$SNP.distances <- 1-adegenet::propShared(x)  
  return(res)
}
