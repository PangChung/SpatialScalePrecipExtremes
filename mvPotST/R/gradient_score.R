#' Gradient score function for the Brown--Resnick model.
#'
#' Compute the peaks-over-threshold gradient score function for the Brown--Resnick model.
#'
#' The function computes the gradient score based on the representation developed by Wadsworth et al. (2014).
#' Margins must have been standardized. The weighting function must be differentiable and verify some properties
#' for consistency, see de Fondeville and Davison (2018) for more details.
#'
#' @author Raphael de Fondeville and modified by Peng Zhong
#' @param obs List of vectors exceeding an R-threshold, see de Fondeville and Davison (2018) for more details.
#' @param loc Matrix of coordinates as given by \code{expand.grid()}.
#' @param vario.fun Semi-variogram function taking a vector of coordinates as input.
#' @param weightFun Function of weights.
#' @param dWeightFun Partial derivative function of \code{weightFun}.
#' @param ... Parameters for \code{weightFun} and \code{dWeightFun}.
#' @param nCores Number of cores used for the computation
#' @param ST logical: whether the model is temporal non-stationary
#' @references de Fondeville, R. and Davison A. (2018). High-dimensional peaks-over-threshold inference. Biometrika, 105(3), 575-592.
#' @references Wadsworth, J. L. and J. A. Tawn (2014). Efficient inference for spatial extreme value
#'  processes associated to log-Gaussian random functions. Biometrika, 101(1), 1-15.
#' @return Evaluation of the gradient score function for the set of observations \code{obs} and semi-variogram \code{vario}.
#' @examples
#' #Define variogram function
#' vario.fun <- function(h){
#'    1 / 2 * norm(h,type = "2")^1.5
#' }
#'
#' #Define locations
#' loc <- expand.grid(1:4, 1:4)
#'
#' #Simulate data
#' obs <- simulPareto(1000, loc, vario.fun)
#'
#' #Evaluate risk functional
#' sums <- sapply(obs, sum)
#'
#' #Define weighting function
#' weightFun <- function(x, u){
#'  x * (1 - exp(-(sum(x / u) - 1)))
#' }
#'
#' #Define partial derivative of weighting function
#' dWeightFun <- function(x, u){
#' (1 - exp(-(sum(x / u) - 1))) + (x / u) * exp( - (sum(x / u) - 1))
#'}
#'
#' #Select exceedances
#' threshold <- quantile(sums, 0.9)
#' exceedances <- obs[sums > threshold]
#'
#' #Evaluate gradient score function
#' scoreEstimation(exceedances, loc, vario, weightFun = weightFun, dWeightFun, u = threshold)
#' @export
scoreEstimation <- function(obs, loc, vario.fun, 
                            weightFun = NULL, dWeightFun = NULL, 
                            nCores = 1L,  
                            ST=FALSE,
                            ... ){
  #Backwards compatibility for versions 0.1.3 and below
  ellipsis <- list(...)
  if("weigthFun" %in% names(ellipsis) && is.null(weightFun)){
    weightFun <- ellipsis[['weigthFun']]
    ellipsis[['weigthFun']] <- NULL
  }
  if("dWeigthFun" %in% names(ellipsis) && is.null(dWeightFun)){
    dWeightFun <- ellipsis[['dWeigthFun']]
    ellipsis[['dWeigthFun']] <- NULL
  }
  if(is.null(weightFun)){
   stop("`weightFun` argument missing, with no default value.")
  }
  if(is.null(dWeightFun)){
    stop("`dWeightFun` argument missing, with no default value.")
  }
  if(!inherits(obs, "list") || length(obs) < 1 || !inherits(obs[[1]], c("numeric","integer"))){
    stop('`obs` must be a list of vectors.')
  }
  if(!inherits(weightFun, "function")){
    stop('`weightFun` must be a function.')
  }
  if(!inherits(dWeightFun, "function")){
    stop('`dWeightFun` must be a function.')
  }
  if(!is.numeric(nCores) || nCores < 1) {
    stop('`nCores` must a positive number of cores.')
  }

  n <- length(obs)
  dim <- nrow(loc)
  if(!ST){
  gamma <- tryCatch({
    dists <- lapply(1:ncol(loc), function(k) {
      outer(loc[,k],loc[,k], "-")
      })
    computeVarMat <- sapply(1:length(dists[[1]]), function(k){
      h <- rep(0,ncol(loc))
      for(j in 1:ncol(loc)){
        h[j] = dists[[j]][k]
      }
      vario.fun(h)
    })
    matrix(computeVarMat, dim, dim)
  }, warning = function(war) {
    war
  }, error = function(err) {
    stop('The semi-variogram is not valid for the locations provided.')
  })
  gammaOrigin <- apply(loc, 1, vario)

  SigmaS <- (outer(gammaOrigin, gammaOrigin, "+") - gamma)
  
  computeScores = function(i){
    obs.i = .subset2(obs,i)
    ind = !is.na(obs.i)
    dim = sum(ind)
    obs.i = obs.i[ind]
    
    sigmaInv <- MASS::ginv(SigmaS[ind,ind])
    sigma <- diag(SigmaS[ind,ind])
    q <- rowSums(sigmaInv)
    
    A <- sigmaInv - q %*% t(q)/sum(q)
    zeroDiagA <- A
    diag(zeroDiagA) <- 0
    mtp <- 2*q / (sum(q)) + 2 + sigmaInv %*% sigma - (q %*% t(q) %*% sigma)/(sum(q))
    
    gradient <- - 1 / 2 * ((A + t(A)) %*% log(obs.i)) *  (1 / obs.i) - 1/2* (1/obs.i) * mtp
    diagHessian <- - 1 / 2 * diag(A + t(A)) * (1/ obs.i^2) + 1 / 2 * ((A + t(A)) %*% log(obs.i)) * (1/obs.i)^2  + 1/2* (1/obs.i)^2 * mtp
    #hack due to intercept in ellipsis function
    weights <- do.call(what = "weightFun", args = c(ellipsis, x = list(obs.i)))
    dWeights  <- do.call(what = "dWeightFun", args = c(ellipsis, x = list(obs.i)))

    sum(2 * (weights * dWeights) * gradient + weights^2 * diagHessian + 1 / 2 * weights^2 * gradient^2)
  }
  }else{
    computeScores = function(i){
      # prepare the data for time i#
      obs.i = .subset2(obs,i)
      ind = !is.na(obs.i)
      dim = sum(ind)
      obs.i = obs.i[ind]
      # compute the variogram matrix
      gamma <- tryCatch({
        dists <- lapply(1:ncol(loc), function(k) {
          outer(loc[ind,k],loc[ind,k], "-")
        })
        computeVarMat <- sapply(1:length(dists[[1]]), function(k){
          h <- rep(0,ncol(loc))
          for(j in 1:ncol(loc)){
            h[j] = dists[[j]][k]
          }
          vario.fun(h,t=i)
        })
        matrix(computeVarMat, dim, dim)
      }, warning = function(war) {
        war
      }, error = function(err) {
        stop('The semi-variogram is not valid for the locations provided.')
      })
      gammaOrigin <- apply(loc[ind,], 1, vario)
      ## compute the sigma matrix and 
      SigmaS <- (outer(gammaOrigin, gammaOrigin, "+") - gamma)
      sigmaInv <- MASS::ginv(SigmaS)
      sigma <- diag(SigmaS)
      q <- rowSums(sigmaInv)
      
      A <- sigmaInv - q %*% t(q)/sum(q)
      zeroDiagA <- A
      diag(zeroDiagA) <- 0
      mtp <- 2*q / (sum(q)) + 2 + sigmaInv %*% sigma - (q %*% t(q) %*% sigma)/(sum(q))
      ## compute the scores
      gradient <- - 1 / 2 * ((A + t(A)) %*% log(obs.i)) *  (1 / obs.i) - 1/2* (1/obs.i) * mtp
      diagHessian <- - 1 / 2 * diag(A + t(A)) * (1/ obs.i^2) + 1 / 2 * ((A + t(A)) %*% log(obs.i)) * (1/obs.i)^2  + 1/2* (1/obs.i)^2 * mtp
      #hack due to intercept in ellipsis function
      weights <- do.call(what = "weightFun", args = c(ellipsis, x = list(obs.i)))
      dWeights  <- do.call(what = "dWeightFun", args = c(ellipsis, x = list(obs.i)))
      
      sum(2 * (weights * dWeights) * gradient + weights^2 * diagHessian + 1 / 2 * weights^2 * gradient^2)
    }
  }
  
  if(nCores > 1){
    scores <- parallel::mclapply(1:n, computeScores,mc.cores = nCores)
  } else {
    scores <- lapply(1:n, computeScores)
  }
  
  return(sum(unlist(scores))/n)
}
