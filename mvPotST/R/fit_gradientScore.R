#' fit the model based on censored likelihood for Brown-Resnick model with gradient score method
#' 
#' @author Peng Zhong
#' @param obs observation list with same marginal distribution, standard pareto distribution
#' @param loc location matrix, e.g., longitude and latitude after transformation
#' @param par parameters for variogram
#' @param vario the variogram function with parameter h,par, and an extra parameter t representing the time position
#' @param u threshold 
#' @param ST logical, default to be FALSE, should the varigoram include a temporal trend
#' @param maxit maximum iteration for the optimization
#' @param nCores number of cores available for parallelization
#' @param weightFun Function of weights.
#' @param dWeightFun Partial derivative function of \code{weightFun}.
#' @references de Fondeville, R. and Davison A. (2018). High-dimensional peaks-over-threshold inference. Biometrika, 105(3), 575-592.
#' @references Wadsworth, J. L. and J. A. Tawn (2014). Efficient inference for spatial extreme value
#'  processes associated to log-Gaussian random functions. Biometrika, 101(1), 1-15.
#' @examples
#' vario <- function(h,par=c(0,log(3)),t=1){ ##return the semi-variogram
#' #reparametrization
#' alpha = 2/(1+exp(-par[1]));lambda1 = par[2]
#' lambda <- exp(lambda1)
#' val=2*(sqrt(sum(h^2))/lambda)^alpha
#' return(val)
#' }
#' #Define locations
#' loc <- expand.grid(1:10, 1:10)
#' 
#' #Simulate data
#' obs <- simulPareto(1000, loc, vario)
#' 
#' #Evaluate risk functional
#' maxima <- sapply(obs, max)
#' thres <- quantile(maxima, 0.9)
#' ncores = detectCores() - 2
#' #Select exceedances
#' exceedances <- obs[maxima > thres]
#' result = fit.gradientScoreBR(obs=exceedances,loc=loc,init=c(0,log(4)),fixed = c(F,F),vario = vario,u = thres,ST = TRUE,nCores = ncores)
#' @export
fit.gradientScoreBR <- function(obs,
                                loc,
                                init,
                                fixed,
                                vario,
                                u,
                                ST=FALSE,
                                maxit = 1000,
                                nCores=1,
                                weightFun=NULL,
                                dWeightFun=NULL,
                                ...){
    
  if(is.matrix(obs)){ #Not converted to list
    obs <- split(obs, row(obs)) #create list
  }
  
  if(is.null(weightFun) | is.null(dWeightFun)){
    weightFun <- function(x, u){
       x * (1 - exp(-(mean(x / u) - 1)))
    }
      #Define partial derivative of weighting function
    dWeightFun <- function(x, u){
      (1 - exp(-(mean(x / u) - 1))) + (x /u/length(x)) * exp( - (mean(x / u) - 1))
    }
  }
  
  fun <- function(par){
    par2 = init
    par2[!fixed] = par 
    if(ST){
      vario.fun <- function(h,t){
        return(vario(h,par2,t))
      }
    }else{
      vario.fun <- function(h){
        return(vario(h,par2))
      } 
    }
    val = scoreEstimation(obs, loc, vario.fun=vario.fun, weightFun = weightFun, dWeightFun=dWeightFun,u = u,nCores = nCores,ST=ST,...)
    return(val)
  }
  init2 = init[!fixed]
  t1 = proc.time()
  result <- optim(init2,fun,control = list(trace=TRUE,maxit=maxit),hessian=FALSE)
  t2 = proc.time() - t1
  init[!fixed] = result$par
  result$par =  init
  result$time = t2 - t1
  return(result)
}
