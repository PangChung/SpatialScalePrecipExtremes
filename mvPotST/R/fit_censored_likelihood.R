#' fit the model based on censored likelihood for Brown-Resnick model
#' 
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
#' 
#' @examples
#' vario <- function(h,par=c(0,3),t=1){ ##return the semi-variogram
#' ## reparametrization
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
#' result = fit.censoredlikelihoodBR(obs=exceedances,loc=loc,init=c(0,4),fixed = c(F,F),vario = vario,u = thres,ST = TRUE,nCores = ncores)
#' @export
fit.censoredlikelihoodBR <- function(obs,
                                     loc,
                                     init,
                                     fixed,
                                     vario,
                                     u,
                                     ST=FALSE,
                                     maxit = 1000,
                                     nCores=1){
  p <- 499
  latticeRule <- genVecQMC(p, (nrow(loc) - 1))
  primeP <- latticeRule$primeP
  vec <- latticeRule$genVec
  
  if(is.matrix(obs)){ #Not converted to list
      obs <- split(obs, row(obs)) #create list
  }
  
  maxima <- unlist(lapply(obs, function(x){max(u,x,na.rm=TRUE)}))
  exceedances <- obs[maxima > u]
  
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
    val = censoredLikelihoodBR(exceedances, loc, vario.fun, rep(u, nrow(loc)),primeP,vec,nCores=nCores,ST=ST)
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
        
                    