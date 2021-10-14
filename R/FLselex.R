# {{{
# selage 
#
#' coverts harvest() and catch.sel() selectivity at age Sa  
#'
#' @param stock Input FLStock object or FLQuant.
#' @param nyears numbers of last years to compute selectivity
#' @param year option to specify year range, overwrites nyears
#' @param maturity if TRUE selage() extract mat_at_age
#' @return FLQuant of Selectivity-at-age Sa  
#' @export
selage <- function (stock, nyears=3,year=NULL,maturity=FALSE){

nyears = min(dims(stock)$maxyear-dims(stock)$minyear+1,nyears)  
  
if(is.null(year)){yr= (dims(stock)$"maxyear"-nyears+1):dims(stock)$"maxyear"} else {
yr = year  
}
if(class(stock)=="FLStock"){
  if(!maturity){flq = harvest(stock)[,ac(yr)]} else {flq =mat(stock)[,ac(yr)]}
} else { flq = stock[,ac(yr)]}

Sa = apply(flq,1,mean,na.rm=T)/max(apply(flq,1,mean,na.rm=T),na.rm=T)

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
#' approximates Age-at-50-selectivity from Sa 
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
#' @param S50 age-at-50-selectivty
#' @param S95 age-at-95-selectivty
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
  if(is.null(Dmin))  Dmin=min(sel[nrow(sa)],0.7)
  
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
  psel = ifelse(selex.dat$age>=max(Smax),psel_c,psel_a)
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
#' fitselex
#'
#' fits selex selectivity function to F_at_age from FLStock 
#'
#' @param Sa selectivity at age Sa = selage(stock)
#' @param S50 init value age-at-50-selectivty
#' @param S95 init value age-at-95-selectivty
#' @param Smax init value age at peak of halfnormal or top of descending slop
#' @param Dcv init value CV of halfnormal determining the steepness of the descending slope
#' @param Dmin init value height of the descending slop at maximum age
#' @param CVlim sets upper CV bound
#' @param ogive fit logistic 
#' @return list with fit and FLquants selex predictions
#' @export 
fitselex <- function(Sa,S50=NULL,S95=NULL,Smax=NULL,Dcv=NULL,Dmin=NULL,CVlim=0.5,logistic=FALSE){
  
  pars =selexpars(Sa=Sa,S50=S50,S95=S95,Smax=Smax,Dcv=Dcv,Dmin=Dmin)
  
 
  
  imp = c(pars)
  imp[4] = CVlim/3
  lower = imp*0.3
  upper = imp*2
  upper[3] = max(as.data.frame(Sa)$age)*2
  lower[3] = min(as.data.frame(Sa)$age[which(as.data.frame(Sa)$data==1)])
  #upper[5] = max(min(1,imp[5]),0.2)
  upper[4] = CVlim
  lower[4] = 0.05
  lower[5]= 0.0001
  upper[5]= 0.75
  if(logistic){
    lower[3]= 100
    upper[3]= 200
    imp[3] = 150
  }
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
#' @param pars 5 selex parameters of class FLPar 
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

varselex = function(pars,stock,step=0.1,amin=NULL,amax=NULL,
                    nyears=3,type=c("crank","shift","dynamic"),return=c("Pars","Sa")){
selpar=pars  
type = type[1]
return = return[1]
if(is.null(amin)) amin = round(selpar[[1]]*0.7,1)
if(type=="crank"){
  if(is.null(amax)) amax = round(selpar[[2]]*0.9,1)
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


#' as.ogive()
#' 
#' Function to set fbar range to F = max(Fa)  
#' @param object selex FLpar or output from fitselex 
#' @param nyears number of end years for reference
#' @return FLStock
#' @export 
as.ogive <- function(object){
  if(class(object)=="list"){
    object$par[3] = 100
    object$par[4] = 0.05
    object$par[5] = 0.05} 
  
    if(class(object)=="FLPar"){
      object[3] = 100
      object[4] = 0.05
      object[5] = 0.05 
    }
    
  return(object)
} #}}}



#' fbar2f()
#' 
#' Function to set fbar range to F = max(Fa)  
#' @param stock class FLStock
#' @param nyears number of end years for reference
#' @param plim set fbar for ages with Selectivy >= plim (default 0.975)
#' @return FLStock
#' @export 
fbar2f <- function(stock,nyears=3,plim=0.975){
sel = selage(stock,nyears=nyears)
age = dims(stock)[["min"]]:dims(stock)[["max"]]
if(plim==1)range(stock)[6:7] = rep(age[which(sel==max(sel))[1]],2)
if(plim<1) range(stock)[6:7] = c(min(age[which(sel>=plim)]),max(age[which(sel>=plim)]))  
 
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
#' Compute equilibrium quantaties across selectivity curves  
#' @param sel selex FLPars() or Sa FLQuants   
#' @param stock stock object of class FLStock 
#' @param sr spawner-recruitment function FLSR
#' @param Fref reference F for which Bref, Cref etc are computed (default=Fsq)  
#' @param nyears number of years for reference conditions   
#' @return FLBRPs object
#' @export
brp.selex = function(sel,stock,sr=NULL,Fref=NULL,nyears=3,plim=0.975){
  obs=TRUE
  object= sel
  stock = allcatch(stock)
  stock = fbar2f(stock,plim=plim)
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
#' @param sel list of selex parameters of class of FLPars()  
#' @param sr optional spawner-recruitment function FLSR
#' @param byears number of backtest years   
#' @param Fref option to input current F value, specify "catch" or  refpts = c("F0","Fmsy","F0.1","Fspr30","Fsq")   
#' @param nyears number of years for referencNULe conditions   
#' @param plim set fbar for ages with Selectivy >= plim (default 0.975)
#' @return FLStocks object
#' @export
selex.backtest = function(sel,stock,sr=NULL,Fref=NULL,byears=10,nyears=3,plim=0.975){
# merge discards to avoid issues in projections
object = sel
if(class(object)=="FLPars") object = par2sa(object,stock)
# set Fbar range 
#stock=fbar2f(stock)
stock=fbar2f(stock,nyears=nyears,plim=plim)

# use geomean sr if sr = NULL (only growth overfishing)
if(is.null(sr)) sr = fmle(as.FLSR(stock,model=geomean),method="BFGS")

# reference selex
Fsq = fabs(stock,nyears=nyears)
fobs = tail(harvest(stock),byears)
Fobs = apply(fobs,2,max)
yrs = dims(fobs)$minyear:dims(fobs)$maxyear
dy = dims(fobs)$minyear-1
Sobs = selage(stock,nyears=nyears)
# prepare stock structure for backtest
stkf = window(stock,end=dy)
stkf = stf(stkf,byears)


if(is.null(Fref)){
  Fref=an(Fobs)
  ctrl_f <- fwdControl(data.frame(year = yrs,
                                  quant = "f",
                                  value = an(Fref)))
  }

if(Fref[1] =="catch"){
  
  ctrl_f <- fwdControl(data.frame(year = yrs,
                                  quant = "catch",
                                  value = an(catch(stock)[,ac(yrs)])))
}

if(Fref[1]%in%c(c("F0","Fmsy","F0.1","Fspr30","Fsq"))){
  refs = refpts(brp(FLBRP(stock,sr)))
  if(Fref=="F0") Fref = refs["virgin","harvest"]+0.001
  if(Fref=="Fmsy") Fref = refs["msy","harvest"]
  if(Fref=="F0.1") Fref = refs["f0.1","harvest"]
  if(Fref=="Fspr30") Fref = refs["spr30","harvest"]
  if(Fref=="Fsq") Fref = Fsq
  ctrl_f <- fwdControl(data.frame(year = yrs,
                          quant = "f",
                          value = an(Fref)))
}


# Recruitment residuals
rec_res = exp(sr@residuals)
# Current selectivity

out = FLStocks(lapply(object,function(x){
  stk = stock
  harvest(stk)[,ac(yrs)][] = x
  stk = fbar2f(stk,plim=plim)
  stk = fwd(stk, control = ctrl_f, sr = sr,residuals=rec_res)
  stk = window(stk,start=dy)
})) 

out = Map(function(x,y){
  x@name = paste(y)
  x
},out,as.list(object@names))
stock@name = "obs"
out = FLStocks(c(FLStocks(obs=stock),out))

return(out)
}

# Add decrete plot option!
#{{{
# selex.fwd() 
#
#' function to forcast nyears under different selex patterns 
#' @param sel list of selex parameters of class of FLPars()  
#' @param stock stock object of class FLStock 
#' @param sr optional spawner-recruitment function FLSR
#' @param fyears number of forecase years   
#' @param Fref option to input current sel refpts = c("F0","Fmsy","F0.1","Fspr30","Fsq")   
#' @param nyears number of years for reference conditions
#' @param plim set fbar for ages with Selectivy >= plim (default 0.975)
#' @param fbar option to not correct for Fbar
#' @return FLStocks object
#' @export
selex.fwd = function(sel,stock,sr=NULL,fyears=50,Fref=NULL,nyears=3,plim=0.975,fbar=FALSE){
  # merge discards to avoid issues in projections
  object = sel
  if(class(object)=="FLPars") object = par2sa(object,stock)
  # set Fbar range 
  #stock=fbar2f(stock)
  if(!fbar) stock=fbar2f(stock,nyears=nyears,plim=plim)
  
  # use geomean sr if sr = NULL (only growth overfishing)
  if(is.null(sr)) sr = fmle(as.FLSR(stock,model=geomean),method="BFGS")
  
  # reference selex
  if(!fbar) Fsq = fabs(stock,nyears=nyears)
  if(fbar) Fsq = mean(tail(fbar(stock),nyears))
  if(is.null(Fref)){ Fref=an(Fsq)} else {
    refs = refpts(brp(FLBRP(stock,sr)))
    if(Fref=="F0") Fref = refs["virgin","harvest"]+0.001
    if(Fref=="Fmsy") Fref = refs["msy","harvest"]
    if(Fref=="F0.1") Fref = refs["f0.1","harvest"]
    if(Fref=="Fspr30") Fref = refs["spr30","harvest"]
    if(Fref=="Fsq") Fref = Fsq
  } 

  dy = range(stock)[["maxyear"]]
  
  yrs = (dy+1):(dy+fyears)
  Sobs = selage(stock,nyears=nyears)
  # prepare stock structure for backtest
  stkf = stf(stock,fyears)
  ctrl_f <- fwdControl(data.frame(year = yrs,
                                  quant = "f",
                                  value = an(Fref)))

  # Current selectivity
  stkobs = stkf
  harvest(stkobs)[,ac(yrs)][] = Sobs
  if(!fbar)stkobs = fbar2f(stkobs,plim=plim)
  stkobs = fwd(stkobs, control = ctrl_f, sr = sr)
  
  
  out = FLStocks(lapply(object,function(x){
    stk = stkf
    harvest(stk)[,ac(yrs)][] = x
    if(!fbar) stk =  fbar2f(stk,plim=plim)
    stk = fwd(stk, control = ctrl_f, sr = sr)
    stk = window(stk,start=dy)
  })) 
 

  out = Map(function(x,y){
    x@name = paste(y)
    x
  },out,as.list(object@names))
  
  stkobs@name = "obs"
  out = FLStocks(c(FLStocks(obs=stkobs),out))
  
  return(out)
}


