#' @rdname GenVectors
#' @encoding UTF-8
#' @export
HaploVectors <- function(x, pop, dist.model = "N", checkdata = TRUE, log.frequencies = FALSE, method = "euclidean", squareroot.dis = TRUE, choices = c(1, 2), analysis = "none", envir, formula, runs = 999, ...){
	res <- HaploDist(x, dist.model)
	res$call <- match.call()
	if (checkdata) {
	  if (is.null(colnames(pop))) {
	    stop("\n Erro in column names of pop data\n")
	  }
	  match.names <- match(rownames(res$individual.per.haplotype), colnames(pop))
	  if (sum(is.na(match.names)) > 0) {
	    print("There are individuals from x data that are not on pop:", quote = FALSE)
	    print(setdiff(rownames(res$individual.per.haplotype), colnames(pop)))
	    stop("\n Individuals not found on pop matrix\n")
	  }
	  pop <- as.matrix(pop[, match.names, drop = FALSE])
	} else{
	  pop <- as.matrix(pop)
	}
	res$haplotype.per.locality <- pop%*%res$individual.per.haplotype
	if(log.frequencies){
	  res$haplotype.per.locality <- log(res$haplotype.per.locality+1)
	}
	res.eigen <- PCPS::pcps(res$haplotype.per.locality, phylodist = res$haplotype.distances, method = method, squareroot = squareroot.dis)
	res.eigen$call <- NULL
	rownames(res.eigen$values) <- sub("pcps", "haplovector", rownames(res.eigen$values))
	colnames(res.eigen$vectors) <- sub("pcps", "haplovector", colnames(res.eigen$vectors))
	colnames(res.eigen$correlations) <- sub("pcps", "haplovector", colnames(res.eigen$correlations))
	res <- c(res, res.eigen)
	res$scores <- summary(res.eigen, choices = choices)$scores$scores.species
	Analysis <- c("none", "adonis", "glm")
	Analysis <- pmatch(analysis, Analysis)
	if (length(Analysis) != 1 | (is.na(Analysis[1]))) {
	  stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
	}
	if(Analysis==2){
	  analysis <- "adonis2.margin"
	}
	FUN <- PCPS::select.pcpsmethod(analysis)
	if(Analysis!=1 & !is.null(FUN)){
	  if(Analysis == 2){
	    test <- PCPS::matrix.p.sig(res$haplotype.per.locality, phylodist = res$haplotype.distances, method.p = method, sqrt.p = squareroot.dis, FUN = FUN, envir = envir, runs = runs, newname = "haplovector", formula = formula, ...)
	  }
	  if(Analysis == 3){
	    choices.analysis <- PCPS::check.formula(formula, colnames(res$vectors))
	    test <- PCPS::pcps.sig(res$haplotype.per.locality, phylodist = res$haplotype.distances, method = method, squareroot = squareroot.dis, choices = choices.analysis, FUN = FUN, envir = envir, runs = runs, newname = "haplovector", formula = formula, ...)
	  }
	  test$call <- NULL
	  test$PCPS.obs <- NULL
	  names(test) <- sub("null.site", "null.turnover", names(test))
	  names(test) <- sub("null.taxa", "null.divergence", names(test))
	  names(test) <- sub("site.shuffle", "turnover", names(test))
	  names(test) <- sub("taxa.shuffle", "divergence", names(test))
	  res <- c(res, test)
	}
	class(res) <- "GenVectors"
	return(res)
}