#' @rdname mutate.names.matrix.p.null
#' @encoding UTF-8
#' @export
matrix.p.sig.mod <- function (comm, phylodist, FUN, runs = 999, parallel = NULL, ...) 
{
  RES <- list(call = match.call())
  res.pcps.null <- PCPS::matrix.p.null(comm, phylodist, runs = runs, calcpcps = FALSE)
  res.pcps.null <- mutate.names.matrix.p.null(res.pcps.null, "pcps", "haplovector")
  RES$P.obs <- res.pcps.null$P.obs
  statistic.obs <- sapply(list(res.pcps.null$P.obs), FUN = FUN, simplify = FALSE, return.model = TRUE, ...)
  RES$model <- statistic.obs[[1]]$mod.obs
  RES$fun <- FUN
  RES$obs.statistic <- statistic.obs[[1]]$statistic.obs
  newClusters <- FALSE
  if (is.numeric(parallel)) {
    parallel <- parallel::makeCluster(parallel, type = "PSOCK")
    newClusters <- TRUE
  }
  if (!inherits(parallel, "cluster")) {
    statistic.null.taxa <- sapply(res.pcps.null$P.null.taxa, FUN = FUN, simplify = FALSE, ...)
    statistic.null.site <- sapply(res.pcps.null$P.null.site, FUN = FUN, simplify = FALSE, ...)
  }
  else {
    statistic.null.taxa <- parallel::parLapply(parallel, res.pcps.null$P.null.taxa, fun = FUN, ...)
    statistic.null.site <- parallel::parLapply(parallel, res.pcps.null$P.null.site, fun = FUN, ...)
  }
  if (newClusters) {
    parallel::stopCluster(parallel)
  }
  RES$statistic.null.site <- do.call("rbind", statistic.null.site)
  RES$statistic.null.taxa <- do.call("rbind", statistic.null.taxa)
  RES$p.site.shuffle <- as.vector(rbind((apply(sweep(RES$statistic.null.site, 2, RES$obs.statistic, ">="), 2, sum)+1)/(runs + 1)))
  RES$p.taxa.shuffle <- as.vector(rbind((apply(sweep(RES$statistic.null.taxa, 2, RES$obs.statistic, ">="), 2, sum)+1)/(runs + 1)))
  return(RES)
}