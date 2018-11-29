#' @title HaploVectors
#' 
#' @description Function to extract haplotypic eigenvectors and perform null model-based tests.
#' 
#' @encoding UTF-8
#' @import PCPS
#' @param x A set of DNA sequences (as an object of class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}}.
#' @param pop A matrix describing the incidence of each individual (columns) in a given locality (rows).
#' @param dist.model A character string used by the function \code{\link{dist.dna}} to specify the evolutionary model to be used to computes pairwise distances from DNA sequences (default dist.model = "N").
#' @param log.frequencies Logical argument (TRUE or FALSE) to specify if transformation of natural logarithms plus one in haplotype per locality data must be applied  (Default log.frequencies = TRUE).
#' @param squareroot.network Logical argument (TRUE or FALSE) to specify if apply square root on haplotype network (Default squareroot.network = TRUE).
#' @param method Dissimilarity index to apply in matrix P, which describes localities by their haplotypic composition, as accepted by vegdist function in vegan package (Default method = "euclidean").
#' @param squareroot.dis Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index in matrix P (Default squareroot.dis = TRUE).
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1, 2)).
#' @param analysis Type of analysis, partial match to "none", "haplotypic" or "haplovector" (Default analysis = "none"). 
#' @param FUN An object of class function to perform the analysis. See Details and Examples.
#' @param envir A matrix with environmental variables for each population, with variables as columns and localities as rows. See Details and Examples.
#' @param choices.analysis Numeric vector to choose the haplovectors used in analysis. See Details and Examples.
#' @param runs Number of permutations for assessing probability of type I error.
#' @param ... Aditional arguments to function \code{\link{matrix.p.sig}} and \code{\link{pcps.sig}}.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{haplotypes}{A list with haplotypes index that identify each observation that share the same haplotype.} 
#' \item{haplotype.distances}{A matrix with pairwise distances between haplotypes.} 
#' \item{haplotype.per.locality}{A matrix with frequency of each haplotype per locality (\bold{\eqn{W}}).} 
#' \item{haplotype.network}{A matrix with haplotype connections described in the network (\bold{\eqn{D_n}}).} 
#' \item{vectors}{Haplotypic eigenvectors (haplovectors).} 
#' \item{values}{Eigenvalues, relative eigenvalues and cumulative relative eigenvalues.} 
#' \item{correlations}{Correlations between haplotypic eigenvectors and haplotypes.} 
#' \item{P}{Matrix of haplotypic composition (\bold{\eqn{P}}).} 
#' \item{scores.haplotypes}{Scores for biplots.} 
#' \item{model}{The observed model returned by FUN.}
#' \item{fun}{The funtion used.}
#' \item{statistic.null.site}{A matrix with null statistic for site shuffle null model.}
#' \item{statistic.null.taxa}{A matrix with null statistic for taxa shuffle null model.}
#' \item{statistic.obs}{Observed statistic, F value to predefined function.}
#' \item{p.site.shuffle}{The p value for the site shuffle null model.}
#' \item{p.taxa.shuffle}{The p value for the taxa shuffle null model.}
#' @details
#' The argument \emph{analysis} specify the type of analysis performed. When \emph{analysis} is equal 
#' "haplotypic" the analysis is performed in matrix of haplotypic composition (using \code{\link{matrix.p.sig}} 
#' function). The argument \emph{formula} can be specified, where the left hand side gives the resemblance data, 
#' right hand side gives the variables. The resemblance data is internally named \emph{p.dist}, thus formula is 
#' an expression of the form \emph{p.dist ~ model}. If \emph{analysis} is equal "haplovector" it is performed 
#' with haplovector (using \code{\link{pcps.sig}} function). In this case, the argument \emph{formula} can be 
#' specified, where the left hand side gives the vectors used, right hand side gives the variables. The vectors 
#' are internally named sequentially \emph{haplovector.1}, \emph{haplovector.2}, \emph{haplovector.3} and so on. 
#' Thus, formula is an expression of the form \emph{haplovector.1 ~ model}. All functions and methods available 
#' in \code{\link{matrix.p.sig}} or \code{\link{pcps.sig}} could be used here, considering the differences 
#' mentioned above.
#' 
#' @seealso \code{\link{HaploNetDist}}
#' @examples 
#' data(sehv)
#' 
#' HaploVectors(sehv$sehv.fas, sehv$sehv.pi, envir = sehv$sehv.envir)
#' 
#' HaploVectors(sehv$sehv.fas, sehv$sehv.pi, analysis = "haplotypic", FUN = FUN.ADONIS, 
#'              envir = sehv$sehv.envir, formula = p.dist~R, 
#'              method.p = "euclidean", sqrt.p = TRUE, runs = 99)
#'              
#' HaploVectors(sehv$sehv.fas, sehv$sehv.pi, analysis = "haplotypic", FUN = FUN.ADONIS2.global, 
#'              envir = sehv$sehv.envir, formula = p.dist~R, 
#'              method.p = "euclidean", sqrt.p = TRUE, runs = 99)
#'              
#' HaploVectors(sehv$sehv.fas, sehv$sehv.pi, analysis = "haplovector", FUN = FUN.GLM, 
#'              envir = sehv$sehv.envir, choices.analysis = 1, 
#'              formula = haplovector.1~R, runs = 99)
#' 
#' @export
HaploVectors <- function(x, pop, dist.model = "N", log.frequencies = TRUE, squareroot.network = TRUE, method = "euclidean", squareroot.dis = TRUE, choices = c(1, 2), analysis = "none", FUN, envir, choices.analysis, runs = 999, ...){
	res <- HaploNetDist(x, pop, dist.model, checkdata = TRUE)
	res$call <- match.call()
	if(squareroot.network){
	  res$haplotype.network <- sqrt(res$haplotype.network)
	}
	if(log.frequencies){
	  res$haplotype.per.locality <- log(res$haplotype.per.locality+1)
	}
	res.eigen <- PCPS::pcps(res$haplotype.per.locality, phylodist = res$haplotype.network, method = method, squareroot = squareroot.dis)
	res.eigen$call <- NULL
	rownames(res.eigen$values) <- sub("pcps", "haplovector", rownames(res.eigen$values))
	colnames(res.eigen$vectors) <- sub("pcps", "haplovector", colnames(res.eigen$vectors))
	colnames(res.eigen$correlations) <- sub("pcps", "haplovector", colnames(res.eigen$correlations))
	res <- c(res, res.eigen)
	res$scores.haplotypes <- summary(res.eigen, choices = choices)$scores$scores.species
	Analysis <- c("none", "haplotypic", "haplovector")
	analysis <- pmatch(analysis, Analysis)
	if (length(analysis) != 1 | (is.na(analysis[1]))) {
	  stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
	}
	if(analysis!=1){
	  if(analysis == 2){
	    test <- matrix.p.sig.mod(res$haplotype.per.locality, phylodist = res$haplotype.network, FUN = FUN, envir = envir, runs = runs, ...)
	  }
	  if(analysis == 3){
	    test <- pcps.sig.mod(res$haplotype.per.locality, phylodist = res$haplotype.network, method = method, squareroot = squareroot.dis, choices = choices.analysis, FUN = FUN, envir = envir, runs = runs, ...)
	  }
	  test$call <- NULL
	  test$PCPS.obs <- NULL
	  res <- c(res, test)
	}
	class(res) <- "HaploVectors"
	return(res)
	
}