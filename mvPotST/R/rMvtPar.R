#' Simulate Pareto random vectors
#'
#' \code{simulPareto} provides \code{n} replicates of the multivariate Pareto distribution
#' associated to log-Gaussian random function with semi-variogram \code{vario}.
#'
#' The algorithm used here is based on the spectral representation of the Brown--Resnick
#' model as described in Dombry et al. (2015). It provides \code{n} replicates conditioned
#' that \code{mean(x) > 1} on the unit Frechet scale.
#'
#' @param n Number of replicates desired.
#' @param loc Matrix of coordinates as given by \code{expand.grid()}.
#' @param vario Semi-variogram function.
#' @param nCores Number of cores used for the computation
#' @param u threshold for the maximum
#' @param shape is the shape parameter in the marginal Pareto distribution
#' @references Dombry, C., Engelke, S. and M. Oesting. Exact simulation of max-stable processes. Biometrika, 103(2), 303-317.
#' @return List of \code{n} random vectors drawn from a multivariate Pareto distribution with semi-variogram \code{vario}.
#' @examples
#' #Define variogram function
#' vario <- function(h){
#'    1 / 2 * norm(h,type = "2")^1.5
#' }
#'
#' #Define locations
#' loc <- expand.grid(1:4, 1:4)
#'
#' #Simulate data
#' obs <- simulPareto(100, loc, vario,u=1)
#' @export


simulPareto <- function(n, loc, vario, nCores = 1,u,shape=1){
  dim <- nrow(loc)
  if(!is.numeric(nCores) || nCores < 1) {
    stop('`nCores`` must a positive number of cores to use for parallel computing.')
  }
  gamma <- tryCatch({
    dists <- lapply(1:ncol(loc), function(i) {
      outer(loc[,i],loc[,i], "-")
    })
    computeVarMat <- unlist(parallel::mclapply(1:length(dists[[1]]), function(i){
      h <- rep(0,ncol(loc))
      for(j in 1:ncol(loc)){
        h[j] = dists[[j]][i]
      }
      vario(h)
    },mc.cores=nCores))
    matrix(computeVarMat, dim, dim)
  }, warning = function(war) {
    war
  }, error = function(err) {
    stop('The semi-variogram provided is not valide for the provided locations.')
  })
  
  d <- nrow(gamma)
  cov.mat <- (outer(gamma[-1,1],gamma[-1,1], "+") - (gamma[-1,-1]))
  chol.mat <- chol(cov.mat)
  
  fun <- function(i){
  proc <- rep_len(0,d)
  mid <- ceiling(d/2);
  quant <- evd::qgpd(0.8,loc=1,scale=1,shape=shape)
  while(mean(proc) < u){
  proc <- rep_len(0,d)  
  k <- sample(1:dim, 1, replace = TRUE)
  proc[-1] <- t(chol.mat) %*% rnorm(d- 1)
  buffer <- exp(proc - proc[k] - gamma[,k])
  proc <- evd::rgpd(1, loc=1, scale=1, shape=shape) * buffer / mean(buffer)
  proc[proc == 0] = .Machine$double.xmin
  }
  return(proc)
  }
  sims <- parallel::mclapply(1:n,fun,mc.cores = nCores)
  return(sims)
}





