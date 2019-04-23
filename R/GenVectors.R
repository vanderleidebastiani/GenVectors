#' @title HaploVectors, SNPVectors and GenVectors
#' 
#' @description Function to extract haplotypic/genetic eigenvectors and perform null model-based tests.
#' 
#' @encoding UTF-8
#' @import PCPS
#' @include HaploVectors.R
#' @include SNPVectors.R
#' @aliases HaploVectors SNPVectors GenVectors print.GenVectors
#' @param x A set of DNA sequences (class "DNAbin" or "haplotype") as used by the function \code{\link{haplotype}} or a set of individual genotypes (an object of class \code{\link{genind}}) as used by the function \code{\link{propShared}}.
#' @param pop A matrix describing the incidence of each individual (columns) in a given locality (rows).
#' @param dist.model A character string used by the function \code{\link{dist.dna}} to specify the evolutionary model to be used to computes pairwise distances from DNA sequences (default dist.model = "N").
#' @param log.frequencies Logical argument (TRUE or FALSE) to specify if transformation of natural logarithms plus one in haplotype per locality data must be applied  (Default log.frequencies = TRUE).
#' @param checkdata Logical argument (TRUE or FALSE) to check if individual sequences in the pop data follow the same order as in the set of DNA sequences (Default checkdata = TRUE).
#' @param method Dissimilarity index to apply in matrix P, which describes localities by their haplotypic/genetic composition, as accepted by vegdist function in vegan package (Default method = "euclidean").
#' @param squareroot.dis Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index in matrix P (Default squareroot.dis = TRUE).
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1, 2)).
#' @param analysis Type of analysis, partial match to "none", "matrix" or "vector" (Default analysis = "none"). 
#' @param analysis.method Type of method to perform the analysis. See Details and Examples (Default analysis.method = "none"). 
#' @param envir A matrix with environmental variables for each population, with variables as columns and localities as rows. See Details and Examples.
#' @param choices.analysis Numeric vector to choose the haplotypic/genetic vectors used in analysis. See Details and Examples.
#' @param runs Number of permutations for assessing probability of type I error.
#' @param ... Aditional arguments to function \code{\link{matrix.p.sig}} and \code{\link{pcps.sig}}.
#' @param distances Matrix containing genetic distances between individuals.
#' @return A list with: \item{call}{Arguments used.} 
#' \item{haplotypes}{A list with haplotypes index that identify each observation that share the same haplotype.} 
#' \item{haplotype.distances}{A matrix with pairwise distances between haplotypes.} 
#' \item{SNP.distances}{A matrix with pairwise distances between alleles.} 
#' \item{individual.per.haplotype}{A matrix with individuals per haplotype.} 
#' \item{haplotype.per.locality}{A matrix with frequency of each haplotype per locality (\bold{\eqn{W}}).} 
#' \item{vectors}{Haplotypic/genetic eigenvectors (haplovectors, snpvectors or geneticvectors).} 
#' \item{values}{Eigenvalues, relative eigenvalues and cumulative relative eigenvalues.} 
#' \item{correlations}{Correlations between haplotypic/genetic eigenvectors and haplotypes/alleles.} 
#' \item{P}{Matrix of haplotypic/genetic composition (\bold{\eqn{P}}).} 
#' \item{scores}{Scores for biplots.} 
#' \item{model}{The observed model returned by FUN.}
#' \item{fun}{The funtion used.}
#' \item{statistic.null.site}{A matrix with null statistic for site shuffle null model.}
#' \item{statistic.null.taxa}{A matrix with null statistic for taxa shuffle null model.}
#' \item{statistic.obs}{Observed statistic, F value to predefined function.}
#' \item{p.site.shuffle}{The p value for the site shuffle null model.}
#' \item{p.taxa.shuffle}{The p value for the taxa shuffle null model.}
#' @details
#' HaploVectors and SNPVectors are two complementary functions to extract haplotypic/genetic eigenvectors 
#' and perform null model-based tests. 
#' 
#' HaploVectors function is based in set of DNA sequences (class "DNAbin" or "haplotype") as used by the 
#' function \code{\link{haplotype}} and pairwise distances from DNA sequences as used in \code{\link{dist.dna}}. 
#' SNPVectors function is based in set of individual genotypes (class "genind") as pairwise distances 
#' between allelesused based in \code{\link{propShared}} function.
#' 
#' The function HaploVectors extract haplotypic eigenvectors and perform null model tests. The argument \emph{analysis} 
#' specify the type of analysis performed. When \emph{analysis} is equal "matrix" the analysis is performed 
#' in matrix of haplotypic composition (using \code{\link{matrix.p.sig}} function). The argument \emph{formula} 
#' can be specified, where the left hand side gives the resemblance data, right hand side gives the variables. 
#' The resemblance data is internally named \emph{p.dist}, thus formula is an expression of the form 
#' \emph{p.dist ~ model}. If \emph{analysis} is equal "vector" it is performed with haplovector 
#' (using \code{\link{pcps.sig}} function). In this case, the argument \emph{formula} can be 
#' specified, where the left hand side gives the vectors used, right hand side gives the variables. The vectors 
#' are internally named sequentially \emph{haplovector.1}, \emph{haplovector.2}, \emph{haplovector.3} and so on. 
#' Thus, formula is an expression of the form \emph{haplovector.1 ~ model}. 
#' 
#' The function SNPVectors work same way, however extract genetic eigenvectors based in distances between alleles. 
#' Similarly the argument \emph{analysis} specify the type of analysis performed. When \emph{analysis} is equal 
#' "matrix" the analysis is performed in matrix of genetic composition and the argument \emph{formula} can be 
#' specified in the same way the HaploVectors function.  If \emph{analysis} is equal "vector" it is performed 
#' with SNPvector and the argument \emph{formula} can be specified. This case the vectors are internally named 
#' sequentially \emph{SNPvector.1}, \emph{SNPvector.2}, \emph{SNPvector.3} and so on. 
#' Thus, formula is an expression of the form \emph{SNPvector.1 ~ model}. 
#' 
#' A third function, called a GenVectors, is also available. In this case, the matrices of distances 
#' between individuals can be supplied directly. This function work same way that other function, but the vectors 
#' are internally named sequentially \emph{geneticvector.1}, \emph{geneticvector.2}, \emph{geneticvector.3} and so on.
#' 
#' All functions and methods available in \code{\link{matrix.p.sig}} or \code{\link{pcps.sig}} could be used here, 
#' considering the differences mentioned above. 
#' 
#' \strong{The argument analysis.method}
#' 
#' The type of analysis performed by this function is specified using the argument \emph{analysis.method}. The current version 
#' of package includes ten predefined function, however additional small functions can be easy specify. All
#' this function uses the environmental variables to analyze the association between haplotypic/genetic
#' composition and environmental predictors. For haplotypic/genetic matrix analysis the predefined 
#' functions available are \emph{"mantel"}, \emph{"adonis"}, \emph{"adonis2.global"} and \emph{"adonis2.margin"}. For haplovectors and SNPvectors
#' analysis the predefined functions available are \emph{"glm"}, \emph{"rda"}, \emph{"gls.marginal"}, 
#' \emph{"gls.sequential"}, \emph{"lme.marginal"} and \emph{"lme.sequential"}.
#' 
#' 
#' @seealso \code{\link{HaploNetDist}}, \code{\link{SNPDist}}, \code{\link{matrix.p.sig}}, \code{\link{pcps.sig}}
#' @examples 
#' data(segv)
#' 
#' HaploVectors(segv$segv.fas, segv$segv.pi, envir = segv$segv.envir)
#' 
#' HaploVectors(segv$segv.fas, segv$segv.pi, analysis = "matrix", 
#'              analysis.method = "adonis", envir = segv$segv.envir, 
#'              formula = p.dist~R, runs = 99)
#' 
#' HaploVectors(segv$segv.fas, segv$segv.pi, analysis = "matrix", 
#'              analysis.method = "adonis2.global", envir = segv$segv.envir, 
#'              formula = p.dist~R, runs = 99)
#' 
#' HaploVectors(segv$segv.fas, segv$segv.pi, analysis = "vector", 
#'              analysis.method = "glm", envir = segv$segv.envir, 
#'              choices.analysis = 1, formula = haplovector.1~R, runs = 99)
#' 
#' 
#' @export
GenVectors <- function(pop, distances, checkdata = TRUE, method = "euclidean", squareroot.dis = TRUE, choices = c(1, 2), analysis = "none", analysis.method = "none", envir, choices.analysis, runs = 999, ...){
  res <- list(call = match.call())
  if (checkdata) {
    if (is.null(rownames(distances))) {
      stop("\n Erro in row names of distances\n")
    }
    if (is.null(colnames(distances))) {
      stop("\n Erro in column names of distances\n")
    }
    if (is.null(colnames(pop))) {
      stop("\n Erro in row names of pop\n")
    }
    match.names <- match(colnames(pop), colnames(distances))
    if (sum(is.na(match.names)) > 0) {
      stop("\n There are individuals from pop data that are not on distances matrix\n")
    }
    distances <- as.matrix(distances[match.names, match.names])
  }
  res$pop <- pop
  res$genetic.distances <- distances
  res.eigen <- PCPS::pcps(res$pop, phylodist = res$genetic.distances, method = method, squareroot = squareroot.dis)
  res.eigen$call <- NULL
  rownames(res.eigen$values) <- sub("pcps", "geneticvector", rownames(res.eigen$values))
  colnames(res.eigen$vectors) <- sub("pcps", "geneticvector", colnames(res.eigen$vectors))
  colnames(res.eigen$correlations) <- sub("pcps", "geneticvector", colnames(res.eigen$correlations))
  res <- c(res, res.eigen)
  res$scores <- summary(res.eigen, choices = choices)$scores$scores.species
  Analysis <- c("none", "matrix", "vector")
  analysis <- pmatch(analysis, Analysis)
  if (length(analysis) != 1 | (is.na(analysis[1]))) {
    stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
  }
  FUN <- PCPS::select.pcpsmethod(analysis.method)
  if(analysis!=1 & !is.null(FUN)){
    if(analysis == 2){
      test <- PCPS::matrix.p.sig(res$pop, phylodist = res$genetic.distances, method.p = method, sqrt.p = squareroot.dis, FUN = FUN, envir = envir, runs = runs, newname = "geneticvector", ...)
    }
    if(analysis == 3){
      test <- PCPS::pcps.sig(res$pop, phylodist = res$genetic.distances, method = method, squareroot = squareroot.dis, choices = choices.analysis, FUN = FUN, envir = envir, runs = runs, newname = "geneticvector", ...)
    }
    test$call <- NULL
    test$PCPS.obs <- NULL
    res <- c(res, test)
  }
  class(res) <- "GenVectors"
  return(res)
}