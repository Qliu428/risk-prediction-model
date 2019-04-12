#########################################################################################################################
### Landmark PSH supermodel sample code              
### Author: Qing Liu
### Reference Paper: Qing Liu, Gong Tang, Joseph P. Costantino, and Chung-Chou H. Chang. "Landmark proportional 
###                  subdistribution hazards models for dynamic prediction of cumulative incidence functions"
#########################################################################################################################


###============== Packages ==================##
library(survival)
require(MASS)
require(dynpred)

###=======================================     Functions   ========================================####

################################################################################################
##### To transform competing risks data into counting process with landmarking ######
################################################################################################
Transdata.LM <- function(tlm,data){
  if(2%in%data$failure){
    d <- data
    
    status <- d$status <- !(d$failure == 0)
    event <- d$failure
    failcode <- 1
    t <- d$time
    
    g <- survfit(Surv(time, status == 0)~1, data=d)        # survival distribution for censoring
    gs <- g$surv
    gt <- g$time
    uct <- sort(unique(t[status==0]))                      # unique censoring time
    uft <- sort(unique(t[event==failcode]))                # unique main event time
    gft <- evalstep(time=gt,stepf=gs,newtime=uft,subst=1)       # numerator of weight          
    
    i.cr <- d$failure == 2
    idcr <- d$id[i.cr]
    tcr <- t[i.cr]                                              # competing risks time 
    gcr <- evalstep(time=gt,stepf=gs,newtime=tcr,subst=1)       # denominator of weight 
    
    # expanding competing risks time
    crdat <- mapply(idcr,tcr,gcr,FUN=function(idcr,tcr,gcr){
      tcr.ft <- uft[tcr<uft]            #failure times larger than cr time 
      tcr.ct <- uct[tcr<=uct]           #censoring times larger than cr time   
      if(length(tcr.ft)>0&length(tcr.ct)>0){
        i.w <- uft>min(tcr.ct)                      #indicator to weight time point
        if(sum(i.w)>0){
          minuftw_1 <- min(which(i.w))-1            #uft prior to first time point
          i.w_1 <- which(i.w) #index of time points to weight
          if(minuftw_1>0) {i.w_1 <- c(minuftw_1,i.w_1)}
          t.w <- uft[i.w_1]
          w <- gft[i.w]/gcr 
          if(minuftw_1>0) {w <- c(1,w) }
          twL <- length(t.w)
          t0 <- c(0,t.w[1:(twL-1)]) 
          t1 <- t.w
        } else {t0<-0; t1<-tcr; w=1}
      } else {t0<-0; t1<-tcr; w=1}
      i.crdat <- data.frame(idcr,t0,t1,w)
      names(i.crdat) <- c("id","start","stop","weight")
      i.crdat },SIMPLIFY=F)
    
    crdat <- do.call("rbind",crdat)  
    newdat <- merge(d,crdat,by="id",all=TRUE) 
    d.weight <- within(newdat,
                       {start[is.na(start)] <- 0
                       stop[is.na(stop)] <- eval(time,newdat)[is.na(stop)]
                       start[start==0] <- tlm
                       weight[is.na(weight)] <- 1})
  }else{
    d.weight <- data
    d.weight$status <-(data$failure!=0)
    d.weight$start <- tlm;d.weight$stop <- data$time; d.weight$weight <- 1}
  return(d.weight)
}

############################################################################################################################
######  To calculate conditional CIF: need to specify the functional forms of s0 for the selected covariates #####
############################################################################################################################
Fwpredict <- function(bet, sf, newdata, tt)
{
  nt <- length(tt)
  sf <- data.frame(time=sf$time,surv=sf$surv,Haz=-log(sf$surv))
  Fw <- data.frame(time=tt,Fw=NA)
  for (i in 1:nt) {
    sfi <- sf
    tti <- tt[i]
    sfi$Haz <- sfi$Haz * exp(newdata[1]*bet[1] + newdata[2]*bet[2]*f1(tt[i]) + newdata[2]*bet[3]*f2(tt[i]) + 
                               newdata[2]*bet[4]*f3(tt[i]) + bet[5]*g1(tt[i]) + bet[6]*g2(tt[i])) ### need to customize this model form according to the landmark PSH supermodel
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w),subst=0)
    Fw$Fw[i] <- 1-exp(-(tmp[2]-tmp[1]))
  }
  return(Fw)
}


###====================================================     Example   ==============================================================####

############################
###### load data  ######
############################
data=read.csv("Example data.csv")

#########################################################################################
###### Specify prediction window and landmark points used to fit supermodel  ######
#########################################################################################
w <- 0.4                                  ## prediction window
s0 <- seq(0,4,by=0.1)                     ## used to fit landmark PSH supermodel
n=length(unique(data$id))                 ## total sample size, number of subjects

###################################################################
###### Build the stacked super dataset and add IPCW ######
###################################################################
## At baseline s0=0
w.LMdata=NULL
LMdata0 <- cutLM(data=data,outcome=list(time="time",status="failure"),
                 LM=s0[1],horizon=s0[1]+w,covs=list(fixed=c("id","z"),
                 varying="zt"),format="long",id="id",rtime="zt.time")
LMdata0$zt.time <- 0
LMdata0$zt <- data$zt[data$zt.time==0]
w.LMdata<- Transdata.LM(s0[1],LMdata0)

for (i in 2:length(s0)){
  LMdata <- cutLM(data=data,outcome=list(time="time",status="failure"),
                  LM=s0[i],horizon=s0[i]+w,covs=list(fixed=c("id","z"),
                  varying="zt"),format="long",id="id",rtime="zt.time")
  w.LMdata.temp <- Transdata.LM(s0[i],LMdata)
  w.LMdata <- rbind(w.LMdata,w.LMdata.temp)
}

#################################################################################################
###### Smoothing on landmark: create parametric functions of landmark time points  ######
#################################################################################################
f1 <- function(t) 1
f2 <- function(t) t
f3 <- function(t) t^2

g1 <- function(t) f2(t)
g2 <- function(t) f3(t)

w.LMdata$zf1 <- w.LMdata$z*f1(w.LMdata$LM)
w.LMdata$zf2 <- w.LMdata$z*f2(w.LMdata$LM)
w.LMdata$zf3 <- w.LMdata$z*f3(w.LMdata$LM)

w.LMdata$ztf1 <- w.LMdata$zt*f1(w.LMdata$LM)
w.LMdata$ztf2 <- w.LMdata$zt*f2(w.LMdata$LM)
w.LMdata$ztf3 <- w.LMdata$zt*f3(w.LMdata$LM)

w.LMdata$LM1 <- g1(w.LMdata$LM)
w.LMdata$LM2 <- g2(w.LMdata$LM)

###### Wald test to assess covariates of which the effects are dependent on the landmark points ####
##test z
LMpsh.z <- coxph(Surv(start,stop,failure==1) ~ zf1 + zf2 + zf3 + zt + LM1 + LM2                 
                + cluster(id), data=w.LMdata, method="breslow",weights=weight)
bet.z <- LMpsh.z$coef
sig.z <- LMpsh.z$var

wh <- 2:3
wald <- t(bet.z[wh]) %*% solve(sig.z[wh,wh]) %*% bet.z[wh]
pval <- 1-pchisq(wald,df=length(wh))
print(data.frame(wald=wald,pval=pval))

##test zt
LMpsh.zt <- coxph(Surv(start,stop,failure==1) ~ ztf1 + ztf2 + ztf3 + z + LM1 + LM2                 
                 + cluster(id), data=w.LMdata, method="breslow",weights=weight)
bet.zt <- LMpsh.zt$coef
sig.zt <- LMpsh.zt$var

wh <- 2:3
wald <- t(bet.zt[wh]) %*% solve(sig.zt[wh,wh]) %*% bet.zt[wh]
pval <- 1-pchisq(wald,df=length(wh))
print(data.frame(wald=wald,pval=pval))


#################################################
###### Fit the landmark PSH supermodel ######
#################################################
LM_PSH_super <- coxph(Surv(start,stop,failure==1) ~ z + ztf1 + ztf2 + ztf3 + LM1 + LM2                 
                + cluster(id), data=w.LMdata, method="breslow",weights=weight)

bet <- LM_PSH_super$coef
sig <- LM_PSH_super$var

#############################################################################################################
###### Predict the conditional CIF Pr{s<T<=s+w, failure=1|T>s, Z(s)} at selected landmark time s ######
#############################################################################################################
s <- c(2.0, 2.2, 2.4)                     ## landmark time points in dynamic prediction, w=0.4

ndata <- data.frame(z=0, ztf1=0, ztf2=0, ztf3=0, LM1=0, LM2=0)
sf2 <- survfit(LM_PSH_super, newdata=ndata)

cond.cif <- matrix(rep(1,times=length(s)*n),ncol=length(s)) #row is the subject, column is for each landmark
for(i in 1:n){
  for(j in 1:length(s)){
    temp.zt.i <- data$zt[data$id==i]
    temp.zt.time.i <- data$zt.time[data$id==i]
    newdata <- c(data$z[data$id==i][1], temp.zt.i[max(which(temp.zt.time.i <= s[j]))]) 
    cond.cif[i,j] <- Fwpredict(bet, sf2, newdata, s[j])[,2]
  }
}

