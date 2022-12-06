## functions to compute consective means
consective.mean <- function(pos,vals,n=30){
  if(pos <= n){
    return(vals[pos])
  }else{
    return(mean(vals[(pos-n+1):pos]))
  }
}

## semi-variogram
vario <- function(h,par=c(0,log(100),0),t=1){ ##return the semi-variogram
  ## reparametrization
  alpha = 2/(1+exp(-par[1]));lambda1 = par[2];lambda2 <- par[3]
  lambda <- exp(lambda1+lambda2*reg.t[t])
  val=(sqrt(sum(h^2))/lambda)^alpha
  return(val)
}

## risk functional
rFun <- function(x,u,xi=est.shape.gpd){
  val = mean((x/u)^xi)^{1/xi}
  return(val)
}

## weight functions that used in the gradient scoring method
weightFun <- function(x,u,xi=est.shape.gpd){
  val = x*(1- exp(1-rFun(x,u,xi)))
  return(val)
}

dWeightFun <- function(x,u,xi=est.shape.gpd){
  r = rFun(x,u,xi)
  val = (1 - exp(-(r - 1))) + (x/u/length(x)) * exp(-(r - 1))*r^(1-xi)*(x/u)^(xi-1)
  return(val)
}

## function to identify the season
getSeason <- function(d){
    y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date(paste0(y,"-3-20"),  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date(paste0(y,"-6-20"),  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date(paste0(y,"-9-20"),  format = "%Y-%m-%d") # Fall Equinox
    
    s.id<- ifelse (d > WS | d <= SE, "Winter",
                   ifelse (d > SE & d <= SS, "Spring",
                           ifelse (d > SS & d <= FE, "Summer", "Fall")))
	return(s.id)
  }

## function to idenfity the year
getYear <- function(d){
	y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    YEND <- as.Date(paste0(y+1,"-1-1"),  format = "%Y-%m-%d")
    if(d > WS && d < YEND) return(y+1) else return(y)
 }

## function to compute the tail-correlation range
solve.h.BR <- function(par,temp,logval=FALSE){
    alpha = par[1]
    lambda = par[2] + par[3]*temp
    if(logval){
      h = lambda + (1/alpha) * (log(2*qnorm(1.95/2)^2))
    }else{
      h=exp(lambda)*(2*qnorm(1.95/2)^2)^(1/alpha)
    }
    return(h)
  }


