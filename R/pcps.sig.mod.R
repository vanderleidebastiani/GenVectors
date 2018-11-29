#' @rdname mutate.names.matrix.p.null
#' @encoding UTF-8
#' @export
pcps.sig.mod <- function (comm, phylodist, method = "bray", squareroot = TRUE, FUN, choices, runs = 999, parallel = NULL, ...) 
{
  RES <- list(call = match.call())
  res.pcps.null <- PCPS::matrix.p.null(comm, phylodist, runs = runs, calcpcps = TRUE, adjpcps = TRUE, choices = choices, method = method, squareroot = squareroot)
  res.pcps.null <- mutate.names.matrix.p.null(res.pcps.null, "pcps", "haplovector")
  RES$PCPS.obs <- res.pcps.null$pcps.obs
  statistic.obs <- sapply(list(res.pcps.null$pcps.obs[, choices, drop = FALSE]), FUN = FUN, simplify = FALSE, return.model = TRUE, ...)
  RES$model <- statistic.obs[[1]]$mod.obs
  RES$fun <- FUN
  RES$obs.statistic <- statistic.obs[[1]]$statistic.obs
  newClusters <- FALSE
  if (is.numeric(parallel)) {
    parallel <- parallel::makeCluster(parallel, type = "PSOCK")
    newClusters <- TRUE
  }
  if (!inherits(parallel, "cluster")) {
    statistic.null.site <- sapply(sapply(res.pcps.null$pcps.null.site, function(x, choices) x[,choices,drop = FALSE], simplify = FALSE, choices = choices), FUN = FUN, simplify = FALSE, ...)
    statistic.null.taxa <- sapply(res.pcps.null$pcps.null.taxa.adj, FUN = FUN, simplify = FALSE, ...)
  }
  else {
    statistic.null.site <- parallel::parLapply(parallel, sapply(res.pcps.null$pcps.null.site, function(x, choices) x[,choices,drop = FALSE], simplify = FALSE, choices = choices), fun = FUN, ...)
    statistic.null.taxa <- parallel::parLapply(parallel, res.pcps.null$pcps.null.taxa.adj, fun = FUN, ...)
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