## prepare the setting: season, region, if bootstrap or not, r risk functional...##
args<-commandArgs(TRUE) 
regions <- c("danube","lena","mississippi","murray","yangtze")
init = c(-1.5,4,0);fixed=c(F,F,F)
bootstrap=FALSE;
region.ind = 1;bootstrap.ind = 1;season.ind=1
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
for (arg in args) eval(parse(text = arg))
## load the packages and the data from marginal fit##
library(mvPotST)
library(parallel)
library(evd)
library(mgcv)
library(lubridate)
load("data/marginal_fit.RData")
source("code/extra_functions.R")
## file where the data should be stored ##
file = paste0("data/fit_pot_ST_",season.ind,"_",regions[region.ind],".Rdata") 
ncores = detectCores()
## choose the r risk functional...##
if(norm.ind==1){
  est.shape.gpd <- est.shape.gpd[region.ind]
}else{
  est.shape.gpd <- 1}

#Define locations
# for(region.ind in c(1,3)){
#  for(season.ind in 1:4){
n.skip = 30
loc.trans = loc.list[[region.ind]]
season.list <- sapply(date.ts,getSeason)
ind.season = sapply(season.list, function(x){x[1]==season[season.ind]})
## load the observations and parameters ## 
obs = split(U.data[[region.ind]],row(U.data[[region.ind]]))
obs = subset(obs,ind.season)
no.obs = sapply(obs,function(x){sum(!is.na(x))})
obs[no.obs>1] = lapply(obs[no.obs>1],evd::qgpd,shape=1,loc=1,scale = 1)
reg.t = tep[[region.ind]][ind.season]
r.obs <- suppressWarnings(unlist(lapply(obs,function(x){if(sum(!is.na(x))!=0){rFun(x[!is.na(x)],u=1,
est.shape.gpd)}else{NA}})))
thres = quantile(r.obs[no.obs > 5],0.95,na.rm=TRUE)
## select the exceedances
ind.exc = no.obs > 5 & r.obs > thres 
stopifnot( sum(ind.exc) > 0 )

reg.t = reg.t[ind.exc]
reg.t = (reg.t - mean(reg.t))/sd(reg.t) 
exceedances <- obs[ind.exc]

if(bootstrap){
  while(!file.exists(file)){Sys.sleep(60)}
  load(file,e<-new.env())
  
  param = e$result$par
  init = param
  rm(e)
  init.seed = as.integer((as.integer(Sys.time())/bootstrap.ind + sample.int(10^5,1))%%10^5)
  set.seed(init.seed)
  boot.index = sample(1:length(exceedances),length(exceedances),replace = TRUE)
  
  exceedances = exceedances[boot.index]
  reg.t = reg.t[boot.index]
  file = paste0("data/fit_pot_ST_bootstrap_",init.seed,"_",season.ind,"_",regions[region.ind],".Rdata")
}

result = fit.gradientScoreBR(obs=exceedances,loc=loc.trans,init=init,fixed = fixed,
vario = vario,u = thres,ST = ST,nCores = ncores,weightFun = weightFun,dWeightFun = dWeightFun)

save(ind.exc,code,result,exceedances,reg.t,est.shape.gpd,thres,file=file)


