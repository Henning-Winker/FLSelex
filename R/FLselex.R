# {{{
# selage 
#
#' coverts harvest() and catch.sel() selectivity at age Sa  
#'
#' @param stock Input FLStock object.
#' @param nyears numbers of last years to compute selectivity
#' @param year option to specify year range, overwrites nyears
#' @return FLQuant of Selectivity-at-age Sa  
#' @export
selage <- function (stock, nyears=3,year=NULL){
if(is.null(year)){yr= (range(stock)["maxyear"]-nyears+1):range(stock)["maxyear"]} else {
yr = year  
}
Sa = apply(harvest(stock[,ac(yr)]),1,mean,na.rm=T)/max(apply(harvest(stock[,ac(yr)]),1,mean,na.rm=T),na.rm=T)
Sa@units = "NA"
return(Sa)
}


#{{{
# fabs() 
#
#' Compute instantaneous F, such that  F_a = Fabs*Sa   
#'
#' @param stock Input FLStock object.
#' @param nys numbers of last years to compute selectivity
#' @return value Fabs  
#' @export
fabs <- function (stock, nyears=3,year=NULL){
  if(is.null(year)){yr= (range(stock)["maxyear"]-nyears+1):range(stock)["maxyear"]} else {
    yr = year  
  }
  Fabs = mean(apply(harvest(stock[,ac(yr)]),2,max,na.rm=T))
  return(Fabs)
}
# }}}


#{{{
# s50
#
#' approximates Age-at-50%-selectivity from Sa 
#'
#' @param Sa selectivity at age Sa = selage(stock)
#' @return s50 
#' @export 
s50 = function(Sa){ 
  Sa = as.data.frame(Sa)
  sao =which(Sa$data>0.5)[1]
  if(sao>1){
    bin = Sa[(sao-1):sao,]
    reg = an(lm(data~age,bin)$coef)
    S50 = max((0.5-reg[1])/reg[2],0.5)
  } else {
    S50 = Sa[1,"data"]/2  
  }
  return(S50)
}


#{{{
# selexpars
#
#' computes initial values for selex pars 
#'
#' @param Sa selectivity at age Sa = selage(stock)
#' @param S50 age-at-50%-selectivty
#' @param S95 age-at-95%-selectivty
#' @param Smax age at peak of halfnormal or top of descending slop
#' @param Dcv CV of halfnormal determining the steepness of the descending slope
#' @param Dmin height of the descending slop at maximum age
#'
#' @return vector of selectivity pars 
#' @export 
selexpars <- function(Sa,S50=NULL,S95=NULL,Smax=NULL,Dcv=NULL,Dmin=NULL){
  S50proxy =s50(Sa)
  sa = as.data.frame(Sa)
  
  age = sa$age
  sel = sa$data
  
  selex.dat = data.frame(age=c((0:3)[which(!0:3%in%age)],age),
                         sel=c(rep(0.01,3)[which(!0:3%in%age)],sel))
  peak = which(selex.dat[,2]>0.85)[1]
  neg=sel[age>peak]
  dcv = ifelse(length(neg)>1, abs(0.5+2*quantile(neg[-1]-neg[-length(neg)],0.2)[[1]]),0.3)
   
  if(is.null(S50))  S50=  S50proxy
  if(is.null(S95))  S95=an(quantile(c(S50proxy,selex.dat[peak[1],1]),0.9))
  if(is.null(Smax))  Smax=selex.dat[max(peak),1]*1.1
  if(is.null(Dcv))  Dcv=0.2 #(1-sel[nrow(sa)])#/length(age[age>peak])
  if(is.null(Dmin))  Dmin=sel[nrow(sa)]
  
  pars = FLPar(S50=S50,S95=S95,Smax=Smax,Dcv=Dcv,Dmin=Dmin)
  return(pars)
}
# }}}



#{{{
# selex
#
#' computes selex curves 
#'
#' @param Sa selectivity at age Sa = selage(stock)
#' @param pars selexpars S50, S95, Smax, Dcv, Dmin 
#' @return FLquants selex predictions
#' @export 
selex <- function(Sa,pars){
  sa = as.data.frame(Sa)
  age = sa$age
  sel = sa$data
  S50 = pars[[1]]
  S95 = pars[[2]]
  Smax =pars[[3]]
  Dcv =pars[[4]]
  Dmin =pars[[5]]
  
  selex.dat = data.frame(age=c((0:3)[which(!0:3%in%age)],age),
                         sel=c(rep(0.01,3)[which(!0:3%in%age)],sel))
  subs = which(selex.dat$age%in%age)
  
  
  psel_a = 1/(1+exp(-log(19)*(selex.dat$age-S50)/(S95-S50)))
  psel_b = dnorm(selex.dat$age,Smax,Dcv*Smax)/max(dnorm(selex.dat$age,Smax,Dcv*Smax))
  psel_c = 1+(Dmin-1)*(psel_b-1)/-1
  psel = ifelse(selex.dat$age<Smax,psel_a,psel_c)
  #psel = ifelse(psel_a<0.95,psel_a,psel_c)
  #psel = psel/max(psel)
  #psel = pmin(psel,0.999)
  
  #resids = log(selex.dat$sel)-log(psel)
  
  #fits=data.frame(age=selex.dat$age[subs],obs=selex.dat$sel[subs],fit=psel[subs],logis=psel_a[subs],halfnorm=psel_b[subs],height=psel_c[subs])
  observed = fitted = logis = hnorm = height = Sa
  fitted[]=matrix(psel[subs])
  logis[]=matrix(psel_a[subs])
  hnorm[]=matrix(psel_b[subs])
  height[]=matrix(psel_c[subs])
  
  pred = FLQuants(observed,fitted,logis,hnorm,height)
  pred@names = c("observed","fitted","logis","hnorm","height")  
  
  
 
  return(pred)
}
# }}}




#{{{
# fitselex
#
#' fits selex selectivity function to F_at_age from FLStock 
#'
#' @param Sa selectivity at age Sa = selage(stock)
#' @param S50 init value age-at-50%-selectivty
#' @param S95 init value age-at-95%-selectivty
#' @param Smax init value age at peak of halfnormal or top of descending slop
#' @param Dcv init value CV of halfnormal determining the steepness of the descending slope
#' @param Dmin init value height of the descending slop at maximum age
#' @return list with fit and FLquants selex predictions
#' @export 
fitselex <- function(Sa,S50=NULL,S95=NULL,Smax=NULL,Dcv=NULL,Dmin=NULL,CVlim=0.5){
  
  pars =selexpars(Sa=Sa,S50=S50,S95=S95,Smax=Smax,Dcv=Dcv,Dmin=Dmin)
  
 
  
  imp = c(pars)
  imp[4] = CVlim/3
  lower = imp*0.3
  upper = imp*2
  #upper[5] = max(min(1,imp[5]),0.2)
  upper[4] = CVlim
  lower[4] = 0.05
  lower[5]= 0.0001

  # Likelihood
  jsel.ll = function(par=imp,data=Sa){
    Sa=data
    flimp = FLPar(S50=par[1],S95=par[2],Smax=par[3],Dcv=par[4],Dmin=par[5])
    pred= selex(Sa=Sa,pars=flimp)
   return(sum(((pred$observed)-(pred$fitted))^2))
  }
  
  
  fit = optim(par=imp, fn = jsel.ll,method="L-BFGS-B",lower=lower,upper=upper, data=Sa, hessian = TRUE)
  fit$par = FLPar(S50=fit$par[1],S95=fit$par[2],Smax=fit$par[3],Dcv=fit$par[4],Dmin=fit$par[5])
  fit$name
  fit$fits = selex(Sa,fit$par)
  return(fit)
}
# }}}



#{{{
#' aopt()
#'
#' Function to compute Aopt, the age where an unfished cohort attains maximum biomass  
#' @param stock class FLStock
#' @return FLQuant with annual spr0y  
#' @export
#' @author Henning Winker
aopt<-function(stock,nyears=3){
  object=stock
  age = dims(object)[["min"]]:dims(object)[["max"]]
  survivors=exp(-apply(m(object),2,cumsum))
  survivors[-1]=survivors[-dim(survivors)[1]]
  survivors[1]=1
  expZ=exp(-m(object[dim(m(object))[1]]))
  if (!is.na(range(object)["plusgroup"]))
    survivors[dim(m(object))[1]]=survivors[dim(m(object))[1]]*(-1.0/(expZ-1.0))
  ba = yearMeans(tail((stock.wt(object)*survivors)[-dims(object)[["max"]],],nyears))
  aopt = age[which(ba==max(ba))[1]]
  # Condition that at aopt fish had spawned 1 or more times on average
  aopt =  max(aopt,(which(yearMeans(tail(object@mat,nyears))>0.5)+1)[1])
  # ensure that aopt <= maxage
  aopt = min(aopt,dims(object)[["max"]]-1)
  
  return(aopt)
}
# }}}


#{{{
# varselex
#
#' function to dynamically vary selex parameters 
#'
#' @param selexpar 5 selex parameters of class FLPar 
#' @param stock optional stock object for tuning of age range  
#' @param step step size of change in one or several pars
#' @param amin start of S50
#' @param amax end of S50, required if stock = NULL
#' @param amax end of S50, required if stock = NULL
#' @param nyears end years for computing aopt() as automated amax limit
#' @param type option of selectivity change "crank","shift" or "dynamic"
#' @param return type of returned object FLPars or FLQuants 
#' @return selex parameters (FLPars) 
#' @export 

varselex = function(selpar,stock,step=0.1,amin=NULL,amax=NULL,
                    nyears=3,type=c("crank","shift","dynamic"),return=c("Pars","Sa")){

type = type[1]
return = return[1]
if(is.null(amin)) amin = round(selpar[[1]]*0.7,1)
if(type=="crank"){
  if(is.null(amax)) amax = round(selpar[[2]]*0.95,1)
}
if(type%in%c("shift","dynamic","selective")){
  if(is.null(amax)) amax = min(max(aopt(stock,nyears),selpar[[1]]),selpar[[1]]+6) 
  }
seqi = seq(amin,amax,step)
diff = seqi-selpar[[1]]
if(type=="crank"){
pars = FLPars(lapply(as.list(diff),function(x){
  out = selpar
  out[1]=out[1]+x
  out
}))  
}  
if(type=="dynamic"){
  pars = FLPars(lapply(as.list(diff),function(x){
    out = selpar
    ds = selpar[3]-selpar[2]
    out[1]= out[1]+x
    out[2] = max(selpar[2],out[1]+0.55)# *1.2
    out[3] = max(out[2]+ds,selpar[3])
    out
  }))  
}  
if(type=="shift"){
  pars = FLPars(lapply(as.list(diff),function(x){
    out = selpar
    out[1:3]=out[1:3]+x
    out
  }))  
}  
pars@names = paste0(seqi)
if(return=="Pars"){ 
  rtn = pars} else {
  rtn = FLQuants(lapply(pars,function(x){
    selex(selage(stock),x)$fitted}))
  
}
return(rtn)
}
# }}}

#' allcatch()
#' 
#' Function to assign all discards to landings for FLBRP refs  
#' @param stock class FLStock
#' @return FLStock
#' @export 
allcatch <- function(stock){
landings.n(stock) = catch.n(stock)  
discards.n(stock)[] = 0.000001  
discards.wt(stock) = stock.wt(stock)
landings(stock) = computeLandings(stock)
discards(stock) = computeDiscards(stock)
return(stock)
}

#' fbar2f()
#' 
#' Function to set fbar range to F = max(Fa)  
#' @param stock class FLStock
#' @param nyears number of end years for reference
#' @return FLStock
#' @export 
fbar2f <- function(stock,nyears=3){
sel = yearMeans(tail(catch.sel(stock),nyears)) 
age = dims(stock)[["min"]]:dims(stock)[["max"]]
range(stock)[6:7] = rep(age[which(sel==max(sel))[1]],2)
return(stock)
}


#' par2sa()
#' 
#' Function to convert selex pars to Sa (FLQuants)  
#' @param pars selexpars 
#' @param object FLStock or FLQuant of Sa 
#' @param nyears number of end years for reference
#' @return FLStock
#' @export 
par2sa <- function(pars,object,nyears=3){
  if(class(object)=="FLStock") object=selage(object,nyears)
  if(class(pars)=="FLPar") pars = FLPars(pars)
  out = FLQuants(lapply(pars,function(x)selex(object,x)$fitted))
  return(out)
}


#{{{
# brp.selex() 
#
#' function to do nyears backtest of selex pattern in FLStocks 
#' @param sel selex FLPars() or Sa FLQuants   
#' @param stock stock object of class FLStock 
#' @param sr spawner-recruitment function FLSR
#' @param Fref reference F for which Bref, Cref etc are computed (default=Fsq)  
#' @param nyears number of years for reference conditions   
#' @return FLBRPs object
#' @export
brp.selex = function(sel,stock,sr=NULL,Fref=NULL,nyears=3){
  obs=TRUE
  object= sel
  stock = allcatch(stock)
  stock = fbar2f(stock)
  if(is.null(Fref)) Fref= fabs(stock,nyears)    
 if(class(object)=="FLPars") object = par2sa(object,stock)
 if(is.null(sr)) sr = fmle(as.FLSR(stock,model=geomean),method="BFGS")
 brps =FLBRPs(lapply(object,function(x){
   stk = stock
   stk@harvest[] = x
   brp = brp(FLBRP(fbar2f(stk),sr))
   brp+FLPar(Fref=Fref) 
  }))
 
 if(obs){
 ref=  brp(FLBRP(fbar2f(stock),sr))
 ref@name = "obs"
 ref=ref+FLPar(Fref=Fref)
 brps = FLBRPs(c(FLBRPs(ref),brps))
 } 
 
 return(brps)
} #}}}








#{{{
# selex.backtest() 
#
#' function to do nyears backtest of selex pattern in FLStocks 
#' @param stock stock object of class FLStock 
#' @param pars list of selex parameters of class of FLPars()  
#' @param sr spawner-recruitment function FLSR
#' @param byears number of backtest years   
#' @param nyears number of years for reference conditions   
#' @param quantity observed "f" or "catch" for fwdControl() 
#' @return FLStocks object
#' @export
selex.backtest = function(sel,stock,byears=10,sr=NULL,nyears=5,quantity=c("f","catch")){
# merge discards to avoid issues in projections

# set Fbar range 
stock=fbar2f(stock)
# use geomean sr if sr = NULL (only growth overfishing)
if(is.null(sr)) sr = fmle(as.FLSR(stock,model=geomean),method = c("BFGS"))


# reference selex
Fsq = fabs(stock,nyears=nyears)
fobs = tail(harvest(stock),byears)
Fobs = apply(fobs,2,max)
Cobs = tail(catch(stock),byears)
sel = yearMeans(tail(catch.sel(stock),nyears)) 
sel = sel/max(sel)
yrs = dims(fobs)$minyear:dims(fobs)$maxyear
dy = dims(fobs)$minyear-1

# prepare stock structure for backtest
sfwd = window(stock,end=dy)
sfwd = stf(refstk,byears)
harvest(refstk)[,ac(yrs)] = fref 




selobs = tail(catch.sel(stock),byears)
selref = selobs
selref[] = sel
fref = Fobs%*%selref
Cobs = tail(harvest(stock),byears)
# prepare stock structure for backtest
refstk = window(stock,end=dy)
refstk = stf(refstk,byears)
harvest(refstk)[,ac(yrs)] = Fref 


}




