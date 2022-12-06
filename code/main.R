## this R code is companying the paper 
## “Are spatial precipitation extremes becoming more intense,
##  wider, or both? An extreme-value perspective”

## load packages ##
library(parallel)
library(fields)
library(mgcv)
library(evd)
library(evgam)
library(mvPotST) ## should be installed before running the dependence model
library(lubridate) ## for handling date

## load the data and extra functions and prepare the dataframe ##
regions <- c("danube", "mississippi")
seasons <- c("Winter","Spring","Summer","Fall")
load("data/marginal_fit.RData")
load("data/tempearature_covariate.RData") 	## stores the temperature data 
											## that are derived from Climate Model Outputs 
											## after kriging and correcting for bias 
source("code/extra_functions.R")

## Exploratory data analysis##
library(ggplot2)
elev.danube <- station.df[[1]]$elev
D = length(elev.danube)
Dt = length(date.ts)
data.danube <- data.frame(precip = 
				drop(as.vector(as.matrix(precip.ts.df[[1]][,-c(1:3)][,order(elev.danube)]))),
                            date=rep(date.ts,D),
                            year=rep(year(date.ts),D),
                            yday=rep(yday(date.ts),D),
                            id = rep(1:D,each=Dt))
elev.miss <- station.df[[3]]$elev
D = length(elev.miss)
Dt = length(date.ts)
data.miss <- data.frame(precip = 
							drop(as.vector(as.matrix(precip.ts.df[[3]][,-c(1:3)][,order(elev.miss)]))),
                            date=rep(date.ts,D),
                            year=rep(year(date.ts),D),
                            yday=rep(yday(date.ts),D),
                            id = rep(1:D,each=Dt))
  data.miss <- data.miss[is.na(data.miss$precip) | (data.miss$precip < 1500 & data.miss$precip >=0),]
  data.danube <- data.danube[is.na(data.danube$precip) | (data.danube$precip < 1500 & data.danube$precip >=0), ]
  data <- list(data.danube,data.miss)
  elev <- list(elev.danube[order(elev.danube)],elev.miss[order(elev.miss)])
  rm(data.miss,data.danube)
  mean.func <- function(x,y=TRUE){
    n <- sum(!is.na(x))
    if(y & n >= 52){
      return(mean(x,na.rm=TRUE))
    }
    if(!y & n >= 10){
      return(mean(x,na.rm=TRUE))
    }
    return(NA)
  }
  lab.regions <- c("Danube","Mississippi")
  ## prepare for the data.frame ##
  library(viridis)
  p1.list <- p2.list <- list()
  for(r in 1:2){
  mean.year <- aggregate(data[[r]]$precip,by=list(data[[r]]$year,data[[r]]$id),FUN=mean.func,y=TRUE)
  mean.day <- aggregate(data[[r]]$precip,by=list(data[[r]]$yday,data[[r]]$id),FUN=mean.func,y=FALSE)  
  mean.day$elev <- rep(elev[[r]],each=366)
  names(mean.day) <- c("day","ID","precip","elev")
  mean.year$elev <- rep(elev[[r]],each=51)
  names(mean.year) <- c("year","ID","precip","elev")
  color_breaks = c(0:10,max(mean.year$precip,na.rm=T))
  main = paste0("Yearly average: ",lab.regions[r])

  p1 <- ggplot(mean.year) + geom_tile(aes(x=year,y=ID,fill=precip)) 
  p1 <- p1 + scale_y_continuous(name="Elevation (m)",breaks=seq(1,length(elev[[r]]),length.out=5),
labels=as.character(elev[[r]][seq(1,length(elev[[r]]),length.out=5)],limits=c(0,length(elev[[r]]))))
  p1 <- p1 + scale_fill_gradientn(name="Precipitation",colours=topo.colors(length(color_breaks),rev = TRUE),
		limits=c(0,max(mean.year$precip,na.rm=T)),
        values=scales::rescale(color_breaks, to = c(0,1), 
        from = c(0,max(mean.year$precip,na.rm=T)))) 
  p1 <- p1 + ggtitle(main) + xlab("Year") + 
        theme(axis.text = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = "transparent", # Needed to add the border
                                      color = "transparent",            # Color of the border
                                      size = 0.5))
  color_breaks = c(0:10,max(mean.day$precip,na.rm=T))
  main = paste0("Daily average: ",lab.regions[r])
  p2 <- ggplot(mean.day) + geom_tile(aes(x=day,y=ID,fill=precip)) 
  p2 <- p2 + scale_y_continuous(name="Elevation (m)",breaks=seq(1,length(elev[[r]]),length.out=5),
		labels=as.character(elev[[r]][seq(1,length(elev[[r]]),length.out=5)],
		limits=c(0,length(elev[[r]]))))
  p2 <- p2 + scale_fill_gradientn(name="Precipitation",colours=topo.colors(length(color_breaks),rev=TRUE),
		limits=c(0,max(mean.day$precip,na.rm=T)),
        values=scales::rescale(color_breaks, to = c(0,1), 
        from = c(0, max(mean.day$precip,na.rm=T)))) 
  p2 <- p2 + ggtitle(main) + xlab("Day") + 
    theme(axis.text = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = "transparent", # Needed to add the border
                                      color = "transparent",            # Color of the border
                                      size = 0.5))
  p1.list[[r]] <- p1; p2.list[[r]] <- p2
  }
  pdf("figures/heatmaps.pdf",width = 6,height = 5,onefile = TRUE)
  for(i in 1:2){
    show(p1.list[[i]])
    show(p2.list[[i]])
  }
  dev.off()  

## prepare the dataframe for marginal fit ##
## the pesudo-uniform scores based on this marginal fit 
## is stored in the list U.data
region.ind = 1 # select region for processing
## prepare the time covariate: year and date
y.ind = year(date.ts);y.ind = y.ind - y.ind[1] + 1
d.ind = yday(date.ts)
alt <- station.df[[region.ind]]$elev/1000 ## elevation of the station
lon <- station.df[[region.ind]]$Y 
lat <- station.df[[region.ind]]$X

## compute the consective temperature averages 
tep.covariate <- tep[[region.ind]]
D = length(lat);Dt = length(date.ts) # dimensions
y = c(as.matrix(precip.ts.df[[region.ind]])[,-c(1:3)])
y.thres <- 10;y = y - y.thres ## remove the values that are below y.thres
## generate the data frame for marginal fitting
data.df <- data.frame(y=y,
                      temp = rep(tep.covariate,times=D),
                      day = rep(d.ind,times=D),
                      year = rep(y.ind,times=D),
                      alt = rep(alt,each=Dt),
                      lon = rep(lon,each=Dt),
                      lat = rep(lat,each=Dt),
                      col=rep(1:D,each=Dt),
                      row=rep(1:Dt,times=D))[!is.na(y) & y>0 & y < 1500,]
data.df = data.df[complete.cases(data.df),] ## select the complete dataframe

## start fitting the marginal model ##
## WARNING! Long time to run ##
message("start Gamma fitting")
formula = y ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10)
results.gam = gam(formula,family=Gamma(link="log"),data=data.df)
est.sig2 <- results.gam$sig2;est.mean <- results.gam$fitted.values
#est.scale <- est.sig2/est.mean;est.shape = est.mean/est.scale
est.shape = 1/est.sig2;est.scale <- est.mean/est.shape;
data.df$est.quantile <- qgamma(0.9,shape = est.shape,scale=est.scale)

message("start binominal fitting")
data.df$y.bin <- as.numeric(data.df$y > data.df$est.quantile)
formula.bin = y.bin ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10)
results.bin <- gam(formula.bin,family = binomial(link="logit"),data=data.df)
est.prob.exceed <- fitted(results.bin) ## fitted exceeding probability

message("start GPD fitting")
library(evgam)
data.df$y.gpd <- data.df$y - data.df$est.quantile
data.df.gpd <- data.df[data.df$y.gpd>0,]
formula.gpd = list(y.gpd ~ log(est.quantile) + temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10),~1)
results.gpd <- evgam(formula.gpd,data=data.df.gpd,family="gpd")
est.scale.gpd = exp(fitted(results.gpd)[,1]);est.shape.gpd = fitted(results.gpd)[1,2]

## transform the data to pseudo uniform scores  ##
est.prob <- pgamma(data.df$y,shape=est.shape,scale=est.scale)
est.prob[data.df$y.bin] <- 1 - est.prob.exceed[data.df$y.bin] + 
est.prob.exceed[data.df$y.bin]*pgpd(data.df.gpd$y.gpd,loc = 0,scale = est.scale.gpd,shape = est.shape.gpd)

U <- U.gpd <- matrix(NA,nrow=Dt,ncol=D)
U[cbind(data.df$row,data.df$col)] <- est.prob/(1+1e-10) ## avoid computational issues


## the pesudo-uniform scores based on this marginal fit 
## is stored in the list U.data
load("data/marginal_fit.RData")

## plot the marginal return level ##
  model.selected <- c("AWI","MIROC","NorESM","AVG") #selected climate models with their averages
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","mississippi")
  lab.regions <- c("Danube","Mississippi")
  y.thres=10
  count = 1
  p.list <- list()
  for(r in 1:2){
    if(r==1){
    load(paste0("data/marginal_model_for_",regions[r],".RData"))
	}else{
    load(paste0("data/marginal_model_for_",regions[r],"1.RData"))
	load(paste0("data/marginal_model_for_",regions[r],"2.RData"))
	load(paste0("data/marginal_model_for_",regions[r],"3.RData"))
	}
    for(s in 1:4){
    data = precip.ts.df[[r*2-1]][,-c(1:3)]
    data[data > 1500] = NA
    loc.ind <- which.max(colMeans(data,na.rm = T))
    ### prepare the data frame to predict ###
    data <- rbind(data.mean.obs.hist.day[[r]][[s]],data.mean.tep.245.day[[r]][[s]],
			data.mean.tep.585.day[[r]][[s]])
    data$type = c(rep("Obs",nrow(data.mean.obs.hist.day[[r]][[s]])),
                  rep("SSP 2-4.5",nrow(data.mean.tep.245.day[[r]][[s]])),
                  rep("SSP 5-8.5",nrow(data.mean.tep.585.day[[r]][[s]])))
    days <- data$date
    d.ind = yday(days)
    y.ind = year(days);y.ind = y.ind - y.ind[1] + 1
    alt <- station.df[[r*2-1]]$elev[loc.ind]/1000
    lon <- station.df[[r*2-1]]$Y[loc.ind]
    lat <- station.df[[r*2-1]]$X[loc.ind]
    main = paste0(lab.regions[r],"; ",model.selected[s],"; ","Station ",loc.ind," (", round(lat,2) ,"N°,",round(lon,2),"E°" ,") ")
    D = length(lat);Dt = length(days)
    data.pred <- data.frame(temp = rep(data$tep,times=D),
                            day = rep(d.ind,times=D),
                            year = rep(y.ind,times=D),
                            alt = rep(alt,each=Dt),
                            lon = rep(lon,each=Dt),
                            lat = rep(lat,each=Dt),
                            col=rep(1:D,each=Dt),
                            row=rep(1:Dt,times=D))
    ind = !is.na(data.pred$temp);data.pred <- data.pred[ind,]
    mean.pred  <- exp(predict.gam(results.gam,newdata = data.pred))
    sig2.pred <- results.gam$sig2
    shape.pred = 1/sig2.pred;scale.pred <- mean.pred/shape.pred;
    quantile.pred <- qgamma(0.9,shape = shape.pred,scale=scale.pred)
    data.pred$est.quantile = quantile.pred
    prob.exceed.pred <- predict.gam(results.bin,newdata=data.pred,type="response")
    gpd.pred <- predict(results.gpd,newdata=data.pred)
    scale.gpd.pred = exp(gpd.pred[,1]);shape.gpd.pred = gpd.pred[1,2]
    return.level = (1-1/(100*365)) 
    prob.gpd <- (return.level - (1-prob.exceed.pred))/prob.exceed.pred
    
    return.value <- y.thres + quantile.pred + qgpd(prob.gpd,loc=0,scale=scale.gpd.pred,shape=shape.gpd.pred)
    by.list = list(season=getSeason(days)[ind],year=sapply(days,getYear)[ind],type=data$type[ind])
    data <-aggregate(return.value,by=by.list,FUN=mean)
	data$season <- as.factor(data$season)
	data$type <- as.factor(data$type)
	data$x[data$type == "Obs" & data$season == "Winter" & data$year == 2016] = NA
    p <- ggplot(data,aes(x=year,y=x,color=season,linetype=type)) + geom_line(alpha=1,size=0.6) 
    p <- p + scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
    p <- p +  xlab("Year") + ylab ("Return level (mm)") 
    p <- p + labs(color="Season",linetype="Group") 
    #p <- p + scale_color_manual(values=hcl.colors(4, "Berlin")) 
    p <- p + ggtitle(main)
    p <- p + theme(axis.text = element_text(size=10), 
                   axis.title.x = element_text(size=14), 
                  axis.title.y = element_text(size=14),
                   plot.title = element_text(hjust = 0.5),
                   panel.border = element_rect(fill = "transparent", # Needed to add the border
                                               color = "black",            # Color of the border
                                               size = 1),
                   panel.background = element_rect(fill = "transparent")) 
    
    #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
    p.list[[count]] <- p;count = count + 1
    print(count)
    }
    }

  pdf("figures/return_level_margins.pdf",width = 6,height = 3,onefile = TRUE)
  for(i in 1:8){
    show(p.list[i])
  }
  dev.off()

## plot the tempearature covariates ##
  model.selected <- c("AWI","MIROC","NorESM","AVG")
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","mississippi")
  lab.regions <- c("Danube","Mississippi")
  count = 1;p.list <- list()
  for(r in 1:2){
    for(s in 1:4){
      data <- rbind(data.mean.obs.hist.day[[r]][[s]],data.mean.tep.245.day[[r]][[s]],
					data.mean.tep.585.day[[r]][[s]])
      data$type = c(rep("Obs",nrow(data.mean.obs.hist.day[[r]][[s]])),
                    rep("SSP 2-4.5",nrow(data.mean.tep.245.day[[r]][[s]])),
                    rep("SSP 5-8.5",nrow(data.mean.tep.585.day[[r]][[s]])))
      p <- ggplot(data,aes(x=date,y=tep,group=type,color=type)) + geom_line(alpha=0.7) 
      p <- p  + xlab("Year") + ylab ("Temperature (°C)")
      p <- p + labs(color='Group') 
      p <- p + scale_color_manual(values=hcl.colors(3, "Berlin")) 
      p <- p + ggtitle(paste0(lab.regions[r],"; ",model.selected[s]))
      p <- p + theme(axis.text = element_text(size=10), 
                     axis.title.x = element_text(size=14), 
                     axis.title.y = element_text(size=14),
                     plot.title = element_text(hjust = 0.5),
                     panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                 color = "black",            # Color of the border
                                                 size = 1),
                     panel.background = element_rect(fill = "transparent"))
      #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
      p.list[[count]] <- p;count = count + 1
    }
  }
  pdf("figures/temperature_covariate.pdf",width = 6,height = 3,onefile = TRUE)
  for(i in 1:8){
    show(p.list[i])
  }
  dev.off()

## dependence fit ##
## call the R script to do the dependence fit ##
## WARNING: Few hours time to run ##
System("Rscript code/Dependence_model_Fit.R \"region.ind=1;season.ind=1;bootstrap=FALSE;norm.ind=1\"")

## plot the dependence estimates ##
load("data/result_pot_10_moving.Rdata",e1<-new.env()) ##for theta = xi
load("data/result_pot_10_moving2.Rdata",e2<-new.env()) ##for theta = 1
library(ggplot2)
library(cowplot)
#load("data/results_pot_moving_temp_full.Rdata")
seasons <- c("Winter","Spring","Summer","Fall")
regions.sub <- c("Danube","Mississippi")
para.char <- c("hat(nu)","hat(lambda)[0]","hat(lambda)[1]")
pic.list <- list()
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

estimates<- data <- data.frame(value=NULL,season=NULL,measure=NULL,variable=NULL,region=NULL)

## combine the estimates with the bootstrap estimates into a single dataframe
for(region.ind in 1:2){
  for(i in 1:3){
  data1 = do.call(rbind,lapply((region.ind*4-4+1):(region.ind*4),function(ind){cbind(e1$estimates.bootstrap.list[[ind]][,i],(ind-1) %% 4 + 1)}))
  data2 = do.call(rbind,lapply((region.ind*4-4+1):(region.ind*4),function(ind){cbind(e2$estimates.bootstrap.list[[ind]][,i],(ind-1) %% 4 + 1)}))
  data = rbind(data,data.frame(value=data1[,1],season=data1[,2],measure=1,variable=i,region=region.ind))
  data = rbind(data,data.frame(value=data2[,1],season=data2[,2],measure=2,variable=i,region=region.ind))
  }
  estimates <- rbind(estimates,data.frame(value=c(e1$estimates.list[[1]]),season=rep(1:4,3),measure=1,variable=rep(1:3,each=4),
  region=1))
  estimates <- rbind(estimates,data.frame(value=c(e1$estimates.list[[2]]),season=rep(1:4,3),measure=1,variable=rep(1:3,each=4),
  region=2))
  estimates <- rbind(estimates,data.frame(value=c(e2$estimates.list[[1]]),season=rep(1:4,3),measure=2,variable=rep(1:3,each=4),
  region=1))
  estimates <- rbind(estimates,data.frame(value=c(e2$estimates.list[[2]]),season=rep(1:4,3),measure=2,variable=rep(1:3,each=4),
  region=2))
}

data$season <- factor(data$season); levels(data$season) = seasons
data$measure <- factor(data$measure) 
data$variable <- factor(data$variable);levels(data$variable) = para.char
data$region = factor(data$region); levels(data$region) = regions.sub

estimates$season <- factor(estimates$season); levels(estimates$season) = seasons
estimates$measure <- factor(estimates$measure) 
estimates$variable <- factor(estimates$variable);levels(estimates$variable) = para.char
estimates$region = factor(estimates$region); levels(estimates$region) = regions.sub
dummy = data.frame(x=0,variable=as.factor(para.char[3]))

pic <- ggplot(data=data,aes(x=season,y=value,color=measure)) + stat_summary(fun.data=f,geom="boxplot",position = position_dodge(1))
pic <- pic + scale_x_discrete(labels = seasons) + scale_color_manual(name=expression(paste(theta,"=")),values=hcl.colors(2, "Zissou 1"),labels=c(expression(hat(xi)),"1")) + coord_flip()
pic <- pic + facet_grid(region~variable,scales="free",shrink = TRUE,labeller=label_parsed) + xlab(NULL) + ylab(NULL)
pic <- pic + theme(axis.text = element_text(size=20), 
          plot.title = element_text(hjust = 0.5),
          strip.text.x=element_text(size=20),
          strip.text.y = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20)) 
pic <- pic + geom_point(data=estimates,aes(x=season,y=value,color=measure),size=2,shape=19,position=position_dodge(1))
pic <- pic + geom_hline(data=dummy,aes(yintercept = 0),size=2,linetype="dashed",color="grey")


pdf(file="figures/boxplot_estimates.pdf",width=18,height=8)
show(pic)
dev.off()


## plot the tail-correlation range ##
  model.selected <- c("AWI","MIROC","NorESM","AVG")
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","mississippi")
  lab.regions <- c("Danube","Mississippi")
  count = 1;p.list <- list()
  for(r in 1:2){
    for(s in 1:4){
        data <- rbind(data.mean.obs.hist[[r]][[s]],data.mean.tep.245[[r]][[s]],data.mean.tep.585[[r]][[s]])
        data$type = c(rep("Obs",nrow(data.mean.obs.hist[[r]][[s]])),
                      rep("SSP 2-4.5",nrow(data.mean.tep.245[[r]][[s]])),
                      rep("SSP 5-8.5",nrow(data.mean.tep.585[[r]][[s]])))
        data$Group = as.factor(paste(data$type,data$season))
        data$type = as.factor(data$type)
        p <- ggplot(data,aes(x=year,y=range)) + geom_line(aes(color=season,group=Group,linetype=type),size=0.6,alpha=1) 
        p <- p + scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
        p <- p + scale_x_continuous(breaks = seq(1965,2100,20),labels = seq(1965,2100,20),expand = c(0, 0), limits = c(1965,2100)) + xlab("Year") + ylab ("Logarithmic range")
        p <- p + labs(linetype="Group",color="Season") + guides(linetype=guide_legend(text=c("Obs","SSP 2-4.5","SSP 5-8.5"))) 
        p <- p + ggtitle(paste0(lab.regions[r],"; ",model.selected[s]))
        p <- p + theme(axis.text = element_text(size=10), 
                       axis.title.x = element_text(size=14), 
                       axis.title.y = element_text(size=14),
                       plot.title = element_text(hjust = 0.5),
                       panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                   color = "black",            # Color of the border
                                                   size = 1),
                       panel.background = element_rect(fill = "transparent"))
        #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
        p.list[[count]] <- p;count = count + 1
    }
  }
  
  pdf("figures/tail_correlation_range.pdf",width = 6,height = 3,onefile = TRUE)
  for(i in 1:length(p.list)){
    show(p.list[i])
  }
  dev.off()
  
## Simulations ##
loc = as.matrix(expand.grid((1:50)/5,(1:50)/5),ncol=2)
solve.h.BR <- function(par,logval=FALSE){
    alpha = par[2]
    lambda <- par[1]
    if(logval){
      h = log(lambda) + (1/alpha) * (log(2*qnorm(1.95/2)^2))
    }else{
    h=lambda * (2*qnorm(1.95/2)^2)^(1/alpha)
    }
  return(h)
}
param_mat <- as.matrix(expand.grid(c(2,5,10),1))
dep_range <- apply(param_mat,1,solve.h.BR)
simu <- list()
for(i in 1:nrow(param_mat)){
	param = param_mat[i,]
	vario <- function(h,par=param){ 
  		## reparametrization
  		alpha = par[2]
  		lambda <- par[1]
  		val=(sqrt(sum(h^2))/lambda)^alpha
  	    return(val)
	}
	set.seed(245254325)
    simu[[i]] <- unlist(simulPareto(1,loc=loc,vario=vario,nCores=4,u=1,shape=1/5))
}

pic <- list()
#MaxLimits = evd::qgpd(0.99,loc=1,scale=1,shape=1/5)
for(i in 1:length(simu)){
 MaxLimits = max(simu[[i]])
 data <- data.frame(x=loc[,1],y=loc[,2],z=pmin(simu[[i]],MaxLimits))
 pic1 <- ggplot(aes(x=x,y=y,fill=z),data=data) + geom_tile()
 pic1 <- pic1 + scale_fill_gradientn(name="Value",colors=hcl.colors(12,"BluYl",rev=TRUE,alpha=1),limits=c(0,MaxLimits))
 pic1 <- pic1 + theme(axis.text = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          plot.title = element_text(hjust = 0.5),aspect.ratio=1) + coord_fixed() 
 pic[[i]] <- pic1 + ggtitle(
  label=bquote(paste("Tail-correlation range:", ~.(round(dep_range[i],1)),";",~lambda==.(param_mat[i,1])))
  )
}
library(cowplot)
pdf("figures/simulation.pdf",width=4.5*3,height=4,onefile=TRUE)
combined_pic <- plot_grid(pic[[1]],pic[[2]],pic[[3]],
nrow=1,rel_widths=c(1,1,1),label_y = expression(alpha==0.5),hjust=-1,align="vh")
show(combined_pic)
dev.off()

