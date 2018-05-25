#' @title HaploNetDist 
#' 
#' @description Function to extract haplotypes, compute pairwise distances between haplotypes, build haplotypic networks and compute the frequency of each haplotype per locality.
#' 
#' @encoding UTF-8
#' @importFrom ape dist.dna
#' @importFrom reshape cast
#' @importFrom pegas haploNet
#' @param x A set of DNA sequences (as an object of class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}}.
#' @param pop A matrix describing the incidence of each individual (columns) in a given locality (rows).
#' @param dist.model A character string used by the function \code{\link{dist.dna}} to specify the evolutionary model to be used to compute pairwise distances from DNA sequences (default dist.model = "N").
#' @param checkdata Logical argument (TRUE or FALSE) to check if individual sequences in the pop data follow the same order as in the set of DNA sequences (Default checkdata = TRUE).
#' @param ... Additional arguments to the function \code{\link{dist.dna}}.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{haplotypes}{A list with haplotypes indices that identify each observation sharing the same haplotype.} 
#' \item{haplotype.distances}{A matrix with pairwise distances between haplotypes.} 
#' \item{haplotype.per.locality}{A matrix with the frequency of each haplotype per locality (\bold{\eqn{W}}).} 
#' \item{haplotype.network}{A matrix with haplotype connections described in the network (\bold{\eqn{D_n}}).}
#' @seealso \code{\link{HaploVectors}}
#' @examples 
#' data(sehv)
#' HaploNetDist(sehv$sehv.fas, sehv$sehv.pi)
#' @export
HaploNetDist <- function(x, pop, dist.model = "N", checkdata = TRUE, ...){
	res <- list(call = match.call())
	pop <- as.matrix(pop)
	n.ind.x <- length(x)
	names.ind.x <- names(x)
	if (checkdata) {
		if (is.null(colnames(pop))) {
			stop("\n Erro in column names of pop data\n")
		}
		match.names <- match(names.ind.x, colnames(pop))
		if (sum(is.na(match.names)) > 0) {
			print("There are individuals from x data that are not on pop:", quote = FALSE)
			print(setdiff(names.ind.x, colnames(pop)))
			stop("\n Individuals not found on pop matrix\n")
		}
		pop <- pop[, match.names, drop = FALSE]
	}
	hap.x <- pegas::haplotype(x)
	ind.hap.x <- attr(hap.x, "index")
	n.hap <- length(ind.hap.x)
	names.hap.x <- paste("haplotype", attr(hap.x, "dimnames")[[1]], sep = ".")
	attr(hap.x, "dimnames")[[1]] <- names.hap.x
	names(ind.hap.x) <- names.hap.x
	hap.dist <- as.matrix(ape::dist.dna(hap.x, model = dist.model, ...))
	ih <- matrix(0, n.ind.x, n.hap, dimnames = list(names.ind.x, names.hap.x))
	for(i in 1:n.hap){
		ih[ind.hap.x[[i]],i]<-1
	}
	ph <- pop %*% ih
	net.x.obs <- as.data.frame(pegas::haploNet(hap.x, d = NULL)[])
	n.net.x <- nrow(net.x.obs)
	col.1 <- data.frame(col.1 = c(1:n.hap)) 
	col.2 <- data.frame(col.2 = c(1:n.hap))
	net.frame.x.full <- data.frame(col.1 = c(1:n.hap, net.x.obs[,1], net.x.obs[,2]), col.2 = c(1:n.hap, net.x.obs[,2], net.x.obs[,1]), value = 1)
	net.x.matrix <- as.matrix(reshape::cast(net.frame.x.full, col.1~col.2, margins = FALSE, fill = 0, add.missing = TRUE))
	rownames(net.x.matrix) <- names.hap.x
	colnames(net.x.matrix) <- names.hap.x
	net.x.matrix.dist <- net.x.matrix*hap.dist
	trunc.hap.x.dist <- ifelse(net.x.matrix.dist==0, max(net.x.obs[,3])+1, net.x.matrix.dist)
	haplonet.dist <- trunc.hap.x.dist*(1-diag(n.hap))
	res$haplotypes <- ind.hap.x
	res$haplotype.distances <- hap.dist
	res$haplotype.per.locality <- ph
	res$haplotype.network <- haplonet.dist
	return(res)
}