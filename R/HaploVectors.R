#' @title HaploVectors
#' 
#' @description Function to extracts haplotypic eigenvectors and performing null model-based tests.
#' 
#' @encoding UTF-8
#' @importFrom PCPS pcps pcps.sig matrix.p.sig
#' @importFrom stats as.formula
#' @param x A set of DNA sequences (as an object of class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}}.
#' @param pop A matrix describing the incidence of each individual (columns) in a given locality (rows).
#' @param dist.model A character string used by the function \code{\link{dist.dna}} to specifying the evolutionary model to be used to computes pairwise distances from DNA sequences (default dist.model = "N").
#' @param log.frequencies Logical argument (TRUE or FALSE) to specify if apply transformation of natural logarithms plus one in haplotype per locality data (Default log.frequencies = TRUE).
#' @param squareroot.network Logical argument (TRUE or FALSE) to specify if use square root in haplotype network (Default squareroot.network = TRUE).
#' @param method Dissimilarity index to aplly in matrix P which describes localities by their haplotypic composition, as accepted by vegdist function in vegan package (Default method = "euclidean").
#' @param squareroot.dis Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index in matrix P (Default squareroot.dis = TRUE).
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1, 2)).
#' @param analysis Type of analysis, partial match to "none", "glm", "adonis" (Default analysis = "none").
#' @param envir A matrix with environmental variables for each population, with variables as columns and locality as rows.
#' @param formula An object of class formula quotation marks used in GLM analysis. See Details.
#' @param AsFactors Encode an environmental variable as factor used in GLM analysis. The sequence is the same that in the environmental data matrix. See Details.
#' @param runs Number of permutations for assessing significance.
#' @param ... Aditional arguments to function \code{\link{matrix.p.sig}} and \code{\link{pcps.sig}}.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{haplotypes}{A list with haplotypes index that identify each observation that share the same haplotype.} 
#' \item{haplotype.distances}{A matrix with pairwise distances between haplotypes.} 
#' \item{haplotype.per.locality}{A matrix with frequency of each haplotype per locality (\bold{\eqn{W}}).} 
#' \item{haplotype.network}{A matrix with haplotype connections described in the network (\bold{\eqn{D_n}}).} 
#' \item{vectors}{The haplotypic eigenvectors (haplovectors).} 
#' \item{values}{The eigenvalues, relative eigenvalues and cumulative relative eigenvalues.} 
#' \item{correlations}{Correlations between haplotypic eigenvectors and haplotypes.} 
#' \item{P}{Matrix of haplotypic composition (\bold{\eqn{P}}).} 
#' \item{scores.haplotypes}{Scores for biplot graphics.} 
#' \item{envir.class}{The class of each variable in environmental data in glm.}
#' \item{formula}{The formula used in glm.} 
#' \item{model}{The model, an object of class glm or adonis.}
#' \item{statistic.obs}{Observed F value.}
#' \item{p.site.shuffle}{The p value for the site shuffle null model.}
#' \item{p.taxa.shuffle}{The p value for the taxa shuffle null model.}
#' @details 
#' The item formula is an expression of the form haplovector.1 ~ envir. The response term must be the 
#' haplovector name, for example haplovector.1, haplovector.2, haplovector.12.
#' 
#' The item AsFactors changes a environmental variable for the class \code{\link{factor}}. The 
#' sequence is the same that in the environmental data matrix, not the order in the formula.
#' Use \code{\link{c}} to combine more that one variable.
#' @seealso \code{\link{HaploNetDist}}
#' @examples 
#' data(sehv)
#' HaploVectors(sehv$sehv.fas, sehv$sehv.pi, envir = sehv$sehv.envir)
#' HaploVectors(sehv$sehv.fas, sehv$sehv.pi, envir = sehv$sehv.envir, 
#'              analysis = "glm", runs = 99, formula = haplovector.1~R)
#' @export
HaploVectors<-function(x, pop, dist.model = "N", log.frequencies = TRUE, squareroot.network = TRUE, method = "euclidean", squareroot.dis = TRUE, choices = c(1, 2), analysis = "none", envir, formula, AsFactors=NULL, runs = 999, ...){
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
	Analysis <- c("none", "glm", "adonis")
	analysis <- pmatch(analysis, Analysis)
	if (length(analysis) != 1 | (is.na(analysis[1]))) {
	  stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
	}
	if(analysis!=1){
	  if(analysis == 2){
	    formula.pcps <- Reduce(paste, deparse(formula))
	    formula.pcps <- stats::as.formula(sub("haplovector", "pcps", formula.pcps))	
	    test <- PCPS::pcps.sig(res$haplotype.per.locality, phylodist = res$haplotype.network, envir = envir, analysis = "glm", method = method, squareroot = squareroot.dis, formula = formula.pcps, AsFactors = AsFactors, runs = runs, ...)
	    test$formula <- formula
	  }
	  if(analysis == 3){
	    test <- PCPS::matrix.p.sig(res$haplotype.per.locality, phylodist = res$haplotype.network, envir = envir, analysis = "adonis", method = method, squareroot = squareroot.dis, runs = runs, ...)
	  }
	  test$call <- NULL
	  res <- c(res, test)
	}
	return(res)
	
}