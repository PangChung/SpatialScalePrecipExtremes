% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rMvParST.R
\name{simu_pareto}
\alias{simu_pareto}
\title{Simulate Pareto random vectors with temporal trend}
\usage{
simu_pareto(n, loc, param, vario, sparse.ind, u, nCores = ncores, shape = 1)
}
\description{
\code{simu_pareto} provide \code{n} replicates of the multivariate Pareto distribution asscociated to
log gaussian random function with semi-variogram \code{vario}, which includes a temporal trend.
}
\details{
The function here is based on the function \code{simulPareto} but with a semi-variogram with a temporal trend

@param n the number of temporal replicates
@param loc matrix of coordinates
@param param parameters used in the semi-variogram function
@param vario semi-variogram function with three parameters h, param, and t
@param nCores number of cores used for parallelization
@param u the threshold
@param shape the shape parameter for marginal Pareto distribution
@param sparse.ind the index of non-missing value to enforce sparsity and ease computation burdens
@return list of \code{n} replicates of random samples
@examples
#Define semi-variogram function
vario <- function(h,par=c(0,log(100),0),t=1){
alpha = 2/(1+exp(-par\link{1}));lambda1 = par\link{2};lambda2 <- par\link{3}
lambda <- exp(lambda1+lambda2\emph{reg.t\link{t})
val=2}(sqrt(sum(h^2))/lambda)^alpha
return(val)
}
loc <- expand.grid(1:4, 1:4)
param = c(0,100,0)
reg.t = 1:10
simu_obs = simu_pareto(10,loc=loc,param=param,vario=vario,u=1,nCores=1)
}
