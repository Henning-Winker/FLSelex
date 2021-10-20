
#' r4col
#' @param n number of colors
#' @param alpha transluscency 
#' @return vector of color codes
#' @export
r4col <- function(n,alpha=1){
  # a subset of rich.colors by Arni Magnusson from the gregmisc package
  # a.k.a. rich.colors.short, but put directly in this function
  # to try to diagnose problem with transparency on one computer
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  rich.vector <- apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha=alpha))
  return(rich.vector)
}


#{{{
# plotselex 
#
#' plots mean selectivity at age Sa across selected years 
#'
#' @param sel selexpars FLPars(s) or Sa FLquants or output from fitselex()
#' @param Sa observed selectivity-at-age (FLQuant) or FLStock 
#' @param obs show observations if TRUE
#' @param compounds option to show selex compounds
#' @param legend.title option for customized legend.title for discrete scenarios only
#' @param max.discrete determines max selectivities after which plot legend becomes continous 
#' 
#' @return FLQuant of Selectivity-at-age Sa  
#' @export
plotselex<- function(sel,Sa=NULL,obs=NULL,compounds=FALSE,colours=NULL,legend.title="S50",max.discrete=8){
  if(is.null(obs) & class(sel)=="FLQuants") obs=FALSE
  if(is.null(obs)) obs=TRUE
  if(class(sel)=="FLQuant") sel = FLQuants(sel)  
  if(class(Sa)=="FLStock") Sa = selage(Sa)
  if(class(sel)=="FLQuants" & is.null(Sa)) Sa = sel[[1]]
  flq = ifelse(class(sel)=="FLQuants",TRUE,FALSE)
  object = sel
  if(is.null(colours)){colf = r4col} else {colf = colours}
  if(class(object)=="list"){
    pars = object$par
    Sa = object$fits$observed
  } else {
    pars=object  
  }
  
if(flq){
    seldat=as.data.frame(pars)
    max.discrete=100
    compounds=FALSE
    p = ggplot(as.data.frame(seldat))+
      geom_line(aes(x=age,y=data,colour=qname))+geom_hline(yintercept = 0.5,linetype="dotted")+
    scale_colour_discrete(legend.title)+
    ylab("Selectivity")+xlab("Age")+
      #scale_x_continuous(breaks = 1:100)+scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
      scale_x_continuous(expand = c(0, 0),breaks = 1:100) + scale_y_continuous(expand = c(0, 0),limits=c(0,1.03),breaks = seq(0, 1, by = 0.25))
    
 } else {

  if(class(pars)=="FLPar"){
    pars = FLPars(pars)
    pars@names = paste0(round(pars[[1]][[1]],2))
  }
  if(is.null(Sa)){
    Sa = FLQuant(c(0.01,0.5,rep(1,ceiling(pars[[1]][[2]]*2)-2)),dimnames=list(age=1:ceiling(pars[[1]][[2]]*2)))  
    obs=FALSE
  }
  # predict
  pdat = FLQuant(0.5,dimnames=list(age=seq(dims(Sa)$min,dims(Sa)$max,0.05))) 
  pred = lapply(pars,function(x){
    selex(pdat,x)
  })
   
  if(length(pred)<max.discrete){
  
  if(compounds==TRUE & length(pred)==1){
    seldat = as.data.frame(pred[[1]][c(3:5,2)]) 
    cols=c(rainbow(3),"black")
  }
  if(compounds==FALSE& length(pred)==1){
    seldat = as.data.frame(pred[[1]][2])  
    cols=c("black")
  }
  if(length(pred)>1){
    seldat = FLQuants(lapply(pred,function(x){
      x[["fitted"]]}))
    compounds=FALSE
  }
  # Plot
  p = ggplot(as.data.frame(seldat))+
    geom_line(aes(x=age,y=data,colour=qname))+geom_hline(yintercept = 0.5,linetype="dotted")
  
  if(length(pred)==1){
    p=p + scale_color_manual("Selex",values=cols)
  } else {
    p = p +scale_colour_discrete(legend.title)
  }
  if(obs & length(pred)==1) p = p+geom_point(data=as.data.frame(Sa),aes(x=age,y=data), fill="white",shape=21,size=2)
  if(obs & length(pred)>1) p = p+geom_line(data=as.data.frame(Sa),aes(x=age,y=data),linetype="dashed")
  
  p = p +ylab("Selectivity")+xlab("Age")+
    #scale_x_continuous(breaks = 1:100)+scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
    scale_x_continuous(expand = c(0, 0),breaks = 1:10)+ scale_y_continuous(expand = c(0, 0),limits=c(0,1.03),breaks = seq(0, 1, by = 0.25))
    
  }
  
  
  if(length(pred)>=max.discrete){
    seldat = FLQuants(lapply(pred,function(x){
      x[["fitted"]]}))
    compounds=FALSE
    dat = as.data.frame(seldat)
    dat$S50 = an(as.character(dat$qname))
    if(obs){
    Sobs = as.data.frame(Sa)
    dat$ao = c(Sobs$age,rep(NA,nrow(dat)-nrow(Sobs)))
    dat$so = c(Sobs$data,rep(NA,nrow(dat)-nrow(Sobs)))
    }
    p = ggplot(data=dat,aes(x=age,y=data,group=S50))+    
    geom_line(aes(color=S50))+
    scale_color_gradientn(colours=rev(colf(20)))+
    ylab("Selectivity")+xlab("Age")+
    scale_x_continuous(expand = c(0, 0),breaks=1:100) + scale_y_continuous(expand = c(0, 0),limits=c(0,1.03))+
    
      theme(legend.key.size = unit(1, 'cm'), #change legend key size
            legend.key.height = unit(1, 'cm'),
            legend.text = element_text(size=7),
            legend.key.width = unit(0.6, 'cm'),
            legend.title=element_text(size=9)
      )
    if(obs) p = p+geom_line(aes(x=ao,y=so),linetype="dashed", na.rm=TRUE)
      
    
    }
 } # end of non flq loop
  return(p)  
}
# }}}


#{{{
# plotselage
#
#' plots mean selectivity at age Sa across selected years 
#'
#' @param stock Input FLStock object.
#' @param nyears numbers of last years to compute selectivity
#' @param year specific years (will overwrite nyears)
#' @return FLQuant of Selectivity-at-age Sa  
#' @export

plotselage<- function(stock,nyears=5,year=NULL){
  if(is.null(year)){yr= (range(stock)["maxyear"]-nyears+1):range(stock)["maxyear"]} else {
    yr = year } 
  Sa = as.data.frame(selage(stock,nyears=nyears,year=year))
  p = ggplot(data=(as.data.frame(catch.sel(stock[,ac(yr)]))),aes(x=age,y=data))+
    geom_line(aes(color = factor(year)))+ theme(legend.title=element_blank())+
    ylab("Selectivity")+xlab("Age")+geom_line(data=Sa,aes(x=age,y=data),size=1)
  return(p)  
}
# }}}



#{{{
# plotFselex 
#
#' plots trade-offs between relative Catch and SSB as a function Selectivity   
#'
#' @param brps output from brp.selex() 
#' @param what type of F c("Fref","Fmsy","F0.1"), ref is by default Fcur
#' @return ggplot   
#' @export

plotFselex = function(brps,what =c("Fref","Fmsy","F0.1")){
  # check if per-recruit
  pr = ifelse(length(params(brps[[1]]))>1,FALSE,TRUE) 
  Obs = ifelse(names(brps)[1]=="obs",TRUE,FALSE)
  what = what[1]
  ref = c("Fref","msy","f0.1")[which(c("Fref","Fmsy","F0.1")%in%what)]  
  
  dat = do.call(rbind,lapply(brps,function(x){
    rps = refpts(x)   
    data.frame(F=an(rps[ref,"harvest"]), Catch = an(rps[ref,"yield"]),BB0 = an(rps[ref,"ssb"]/rps["virgin","ssb"])) 
  }))    
  d.=data.frame(sel=brps@names,dat)
  rownames(d.) = 1:nrow(d.)
  d.[d.[]<0] = 0
  dat = d.
  if(Obs)  dat=d.[-1,]
  
  scale = mean(dat$BB0)/min(dat$Catch/max(dat$Catch))
  p <- ggplot(dat, aes(x = an(sel)))
  if(!pr){
  p <- p + geom_line(aes(y = BB0/scale, colour = "SSB"))
  p <- p + geom_line(aes(y = Catch/max(Catch), colour = "Yield"))
  p <- p + scale_y_continuous(sec.axis = sec_axis(~.*scale , name = expression(SSB/SSB[0])))
  } else {
    p <- p + geom_line(aes(y = BB0/scale, colour = "SPR"))
    p <- p + geom_line(aes(y = Catch/max(Catch), colour = "YPR"))
    p <- p + scale_y_continuous(sec.axis = sec_axis(~.*scale , name = expression(SPR/SPR[0])))
  }
  p <- p + scale_colour_manual(values = c("red", "blue"))
  if(!pr) p <- p + labs(y = "Relative Yield",x = "Age-at-50%-Selectivity",colour = "")
  if(pr) p <- p + labs(y = "Relative YPR",x = "Age-at-50%-Selectivity",colour = "")
  p <- p + theme(legend.position = "bottom")
  Aopt = dat[dat$Catch==max(dat$Catch),]
  p = p+geom_segment(x=an(Aopt$sel),xend=an(Aopt$sel),y=0,yend=Aopt$BB0/scale,col="red",linetype="dotted")
  p = p+geom_segment(x=an(Aopt$sel),xend=max(an(dat$sel))+10,y=Aopt$BB0/scale,yend=Aopt$BB0/scale,col="red",linetype="dotted")
  p = p+geom_segment(x=0,xend=an(Aopt$sel),y=1,yend=1,col="blue",linetype="dotted")
  #p = p+geom_segment(x=an(Aopt$sel),xend=an(Aopt$sel),y=Aopt$BB0/scale,yend=1,col="blue",linetype="dotted")
  p=p + geom_point(aes(x= an(Aopt$sel),y=Aopt$Catch/max(Catch)),col="blue",size=3) 
  p=p + geom_point(aes(x= an(Aopt$sel),y=Aopt$BB0/scale),col="red",size=3) 
  
  
  #p = p+annotate("text",x=mean(an(dat$sel)),y=1,label=paste0(what,"=",round(dat$F[1],3)))
  if(Obs){
    S50 = s50(brps[[1]]@landings.sel/max(brps[[1]]@landings.sel))
    p = p+geom_segment(x =S50,xend=S50,y=0,yend=0.95 ,linetype="dashed")
    if(pr & what =="Fmsy") what="Fmax" 
    p = p+annotate("text",x=S50+0.03,y=0.97,label=paste0(what," = ",round(dat$F[1],3)),size=3)
  }
  return(p)
}

#{{{
#' Fselex 
#'
#' Returns Table with trade-offs between relative Catch and SSB as a function Selectivity   
#'
#' @param brps output from brp.selex() 
#' @param what type of F c("Fref","Fmsy","F0.1"), ref is by default Fcur
#' @return data.frame(F,rel.yield,rel.ssb)   
#' @export
Fselex = function(brps,what =c("Fref","Fmsy","F0.1")){
  Obs = ifelse(names(brps)[1]=="obs",TRUE,FALSE)
  what = what[1]
  ref = c("Fref","msy","f0.1")[which(c("Fref","Fmsy","F0.1")%in%what)]  
  dat = do.call(rbind,lapply(brps,function(x){
    rps = refpts(x)   
    data.frame(F=an(rps[ref,"harvest"]), rel.yield = an(rps[ref,"yield"]),rel.ssb = an(rps[ref,"ssb"]/rps["virgin","ssb"])) 
  }))    
  d.=data.frame(sel=brps@names,dat)
  d.$rel.yield= d.$rel.yield/max(d.$rel.yield)
  rownames(d.) = 1:nrow(d.)
  d.[d.[]<0] = 0
  
  #if(Obs)  dat=d.[-1,]
  
  return(d.)
}




#{{{
#' ploteqselex()
#
#' Isopleth plot of Selectivity vs F
#'
#' @param brps output from brp.selex() 
#' @param Fmax upper possible limit of  F range
#' @param panels choice of plots 1:4
#' \itemize{
#'   \item 1  F over Relative Yield
#'   \item 2  F over SRP or SSB (requires ssr)
#'   \item 3  Isopleth Yield
#'   \item 4  Isopleth SPR or SSB (requires ssr)
#' }    
#' @param ncol number of columns
#' @param colours optional, e.g. terrain.col, rainbow, etc.
#' @param Ftrg  option to dynamic Ftrg=c("none","fmsy","f0.1")
#' \itemize{
#'   \item "none"  
#'   \item "msy"  
#'   \item "f0.1"
#' }
#' @return ggplot   
#' @export
ploteqselex = function(brps,Fmax=2.,panels=NULL, ncol=NULL,colours=NULL,Ftrg=c("none","msy","f0.1")){
# Colour function
if(is.null(colours)){colf = r4col} else {colf = colours}
if(is.null(panels)) panels=1:4
if(is.null(ncol)){
  if(length(panels)%in%c(1,3)){ncol=length(panels)} else {ncol=2}
}

# Check range
if(paste(brps[[1]]@model)[3]%in%c("a + ssb/ssb - 1")){
pr = TRUE
lim = min(Fmax,max(2*refpts(brps[[1]])["f0.1","harvest"],refpts(brps[[1]])["Fref","harvest"]*1.5,dims(brps[[1]])[["min"]]))
quants = c("YPR","SPR")
labs = c("Relative YPR",expression(SPR/SPR[0]))

} else {
pr = FALSE
lim = min(Fmax,max(2*refpts(brps[[3]])["msy","harvest"],refpts(brps[[1]])["Fref","harvest"]*1.05,dims(brps[[1]])[["min"]]))
quants = c("Yield","SSB")
labs = c("Relative Yield",expression(SBB/SSB[0]))
}

# Prep some data for plotting
fbar(brps[[1]]) = seq(0,lim,lim/101)[1:101]
obs = data.frame(obs="obs",model.frame(metrics(brps[[1]],list(ssb=ssb, harvest=fbar, rec=rec, yield=landings)),drop=FALSE))
obs[,8:11][obs[,8:11]<0] <- 0
S50 = as.list(an(brps[-1]@names))
isodat = do.call(rbind,Map(function(x,y){
  fbar(x) = seq(0,lim,lim/101)[1:101]
  mf =  model.frame(metrics(x,list(ssb=ssb, harvest=fbar, rec=rec, yield=landings)),drop=FALSE)
  data.frame(S50=y,as.data.frame(mf))  
},brps[-1],S50))
isodat$yield = isodat$yield#/max(isodat$yield)
isodat[,8:11][isodat[,8:11]<0] <- 0
Fobs = an(refpts(brps[[1]])["Fref","harvest"])
Yobs = an(refpts(brps[[1]])["Fref","yield"])
Sobs = an(refpts(brps[[1]])["Fref","ssb"])#/an(refpts(brps[[1]])["virgin","ssb"])
Sa=brps[[1]]@landings.sel/max(brps[[1]]@landings.sel)
S50obs = s50(Sa)
isodat$Fo = c(obs$harvest,rep(NA,nrow(isodat)-nrow(obs)))
isodat$Yo = c(obs$yield,rep(NA,nrow(isodat)-nrow(obs)))/max(isodat$yield)
isodat$So = c(obs$ssb,rep(NA,nrow(isodat)-nrow(obs)))/max(isodat$ssb)

if(Ftrg[1]!="none"){
ftrg = do.call(rbind,Map(function(x,y){
frp= c(refpts(x)[Ftrg,"harvest"])
data.frame(s50=y,ftrg=frp)
},brps[-1],S50))
ftrg = ftrg[ftrg$ftrg<lim,]
}
  # F vs Yield
  P1 = ggplot(data=isodat,aes(x=harvest,y=yield/max(yield),group=S50))+
  geom_line(aes(color=S50))+geom_line(aes(x=Fo,y=Yo),size=0.7,linetype="dashed", na.rm=TRUE)+
  scale_color_gradientn(colours=rev(colf(20)))+ylab(labs[1])+
  geom_segment(aes(x = Fobs, xend = Fobs, y = 0, yend = Yobs/max(yield)), colour = "black",size=0.3,linetype="dotted")+
  geom_segment(aes(x = 0, xend = Fobs, y = Yobs/max(yield), yend = Yobs/max(yield)), colour = "black",size=0.3,linetype="dotted")+
  geom_point(aes(x=Fobs,y=Yobs/max(yield)),size=2)+
  xlab("Fishing Mortality")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'),
          legend.text = element_text(size=7),
          legend.key.width = unit(0.6, 'cm'),
          axis.title=element_text(size=10),
          legend.title=element_text(size=9)
    )+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0),limits=c(0,1))
  # F vs SSB  
  P2 = ggplot(data=isodat,aes(x=harvest,y=ssb/max(ssb),group=S50))+
  geom_line(aes(color=S50))+geom_line(aes(x=Fo,y=So),size=0.7,linetype="dashed", na.rm=TRUE)+
  scale_color_gradientn(colours=rev(colf(20)))+
  geom_segment(aes(x = Fobs, xend = Fobs, y = 0, yend = Sobs/max(ssb)), colour = "black",size=0.3,linetype="dotted")+
  geom_segment(aes(x = 0, xend = Fobs, y = Sobs/max(ssb), yend = Sobs/max(ssb)), colour = "black",size=0.3,linetype="dotted")+
  geom_point(aes(x=Fobs,y=Sobs/max(ssb)),size=2)+
  ylab(labs[2])+xlab("Fishing Mortality")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'),
          legend.key.width = unit(0.6, 'cm'),
          legend.text = element_text(size=7),
          axis.title=element_text(size=10),
          legend.title=element_text(size=9)
    )+ 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0),limits=c(0,1))
  
  # Isopleth plot Yield
  colbr = c(seq(0,0.6,0.2),seq(0.7,0.9,0.1),0.95,1)
  conbr = c(0,0.2,0.4,seq(0.5,0.9,0.1),0.95,1)
  nbr = length(colbr)
  P3=ggplot(isodat, aes(x=harvest,y=S50))+
  geom_raster(aes(fill = yield/max(yield)), 
              interpolate = T, hjust = 0.5, vjust = 0.5)+ 
  metR::geom_contour2(aes(z=yield/max(yield)),color = grey(0.4,1),breaks =conbr )+ 
  metR::geom_text_contour(aes(z=yield/max(yield)),stroke = 0.2,size=3,skip=0,breaks = conbr)+
  scale_fill_gradientn(colours=rev(colf(nbr+3))[-c(10:11,13)],limits=c(-0.03,1), breaks=colbr, name=paste(quants[1]))+
  geom_point(aes(x=Fobs,y=S50obs),size=2.5)+
  geom_segment(aes(x = Fobs, xend = Fobs, y = min(S50), yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
  geom_segment(aes(x = 0, xend = Fobs, y = S50obs, yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'),
          legend.key.width = unit(0.6, 'cm'),
          legend.text = element_text(size=7),
          axis.title=element_text(size=10),
          legend.title=element_text(size=9)
    )+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(-0.03, 0))+
  ylab("Age-at-50%-Selectivity")+xlab("Fishing Mortality")

  # Isopleth SSB
  colbr = c(0.05,seq(0,1,0.1))
  conbr = c(seq(0.05,0.4,0.05),seq(0.5,0.6,1),1)
  nbr = length(colbr)
  P4 = ggplot(isodat, aes(x=harvest,y=S50))+
  geom_raster(aes(fill = ssb/max(ssb)), 
              interpolate = T, hjust = 0.5, vjust = 0.5)+ 
  metR::geom_contour2(aes(z=ssb/max(ssb)),color = grey(0.4,1),breaks =conbr )+
  metR::geom_text_contour(aes(z=ssb/max(ssb)),stroke = 0.2,size=3,skip=0,breaks = conbr)+
  scale_fill_gradientn(colours=rev(colf(nbr+4))[-c(4:7)],limits=c(-0.03,1), breaks=colbr, name=quants[2])+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        legend.text = element_text(size=7),
        axis.title=element_text(size=10),
        legend.title=element_text(size=9)
        )+
        geom_point(aes(x=Fobs,y=S50obs),size=2.5)+
  geom_segment(aes(x = Fobs, xend = Fobs, y = min(S50), yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
  geom_segment(aes(x = 0, xend = Fobs, y = S50obs, yend = S50obs), colour = "black",size=0.3,linetype="dotted")+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(-0.03, 0))+
  ylab("Age-at-50%-Selectivity")+xlab("Fishing Mortality")
  if(Ftrg[1]!="none"){
    P3 = P3+geom_point(data=ftrg,aes(x=ftrg,y=s50),shape = 21, colour = "black",size=1.1, fill = "white")
    P4= P4+geom_point(data=ftrg,aes(x=ftrg,y=s50),shape = 21, colour = "black",size=1.1, fill = "white")
  }
  plots <- list(P1=P1,P2=P2,P3=P3,P4=P4)
  
  if(length(panels)>1) return(gridExtra::grid.arrange(grobs =  plots[panels], ncol = ncol))  
  if(length(panels)==1) return(plots[[panels]])  
} #}}}

#{{{
#' plotprjselex()
#
#' Projection plot of FLStocks
#'
#' @param object FLStocks output from selex.forward/selex.backtest() 
#' @param panels choice of plots 1:10
#' \itemize{
#'   \item 1  spawning biomass (SSB)
#'   \item 2  catch
#'   \item 3  Fjuv/F 
#'   \item 4  Percentage juveniles in catch in numbers
#'   \item 5  Harvest Rate = Catch/Vuln.Bio
#'   \item 6  vuln.bio 
#'   \item 7  fishing mortality F
#'   \item 8  total biomass
#'   \item 9  Frec/F (Vasilakopoulos et al. 2020)
#'   \item 10 Percentage in catch >= aopt
#' }   
#' @param ncol number of columns
#' @param nyears number of years back from assessment year shown
#' @param colours optional, e.g. terrain.col, rainbow, etc.
#' @param max.discrete determines max selectivities after which plot legend becomes continous 
#' @return ggplot   
#' @export

plotprjselex = function(object,panels=NULL, ncol=NULL,colours=NULL,nyears=NULL,max.discrete=10){
  
  discrete = ifelse(length(object)>max.discrete,FALSE,TRUE)
    
  if(is.null(nyears)) nyears = dim(object[[1]])[2]-dim(object[[2]])[2]+1
  if(is.null(colours)){colf = r4col} else {colf = colours}
  if(is.null(panels)) panels=1:6
  if(is.null(ncol)){
    if(length(panels)%in%c(1,3)){ncol=1} else {ncol=2}
  }
  object = window(object,start=min(do.call(c,lapply(object,function(x){dims(object[[1]])[["minyear"]]}))),
           end=max(do.call(c,lapply(object,function(x){dims(object[[1]])[["maxyear"]]}))))
  nc = length(object)+3
  #Prep data
  
  
  
  dat = do.call(rbind,lapply(object,function(x){
  df = as.data.frame(FLQuants(
          Catch = catch(x),
          SSB = ssb(x),
          Fjuv.F = apply(harvest(x)*(1-mat(x)),2,mean)/apply(harvest(x),2,max),
          Prop.Juveniles = apply(catch.n(x)*(1-mat(x)),2,sum)/apply(catch.n(x),2,sum)*100,
          Harvest=catch(x)/apply(stock.wt(x)*stock.n(x)*catch.sel(x),2,sum),
          Vuln.Bio = apply(stock.wt(x)*stock.n(x)*catch.sel(x),2,sum),
          F = fbar(x),
          Biomass = stock(x),
          Frec.F = harvest(x)[1,]/apply(harvest(x),2,max),
          Prop.Aopt  = apply(0.001+catch.n(x)[ac(aopt(x):range(x)[["max"]]),],2,sum)/apply(0.001+catch.n(x),2,sum)*100
          
          )[panels])
  data.frame(sel=x@name,df)
  }))
  
  
  dat$qname = ifelse(dat$qname=="Harvest","Harvest rate",paste(dat$qname))
  dat$qname = ifelse(dat$qname=="Prop.Juveniles","Juvenile catch (%)",paste(dat$qname))
  dat$qname = ifelse(dat$qname=="Frec.F","Frec/F",paste(dat$qname))
  dat$qname = ifelse(dat$qname=="Fjuv.F","Fjuv/F",paste(dat$qname))
  dat$qname = ifelse(dat$qname=="Prop.Aopt","Proportion Aopt (%)",paste(dat$qname))
  
  d. = dat[dat$sel!="obs",]
  obs = dat[dat$sel=="obs",]
  d.$obs = c(obs$data,rep(NA,nrow(d.)-nrow(obs)))
  d.$qname = factor(d.$qname,levels=unique(d.$qname))
  
  if(!discrete){
    d.$S50 = an(d.$sel)
  p=ggplot(d.,aes(x=year,y=data,group=sel))+geom_line(aes(col=S50), na.rm=TRUE)+
  scale_color_gradientn(colours=rev(colf(nc)[-c(1:2)]))+
  facet_wrap(~qname, scales="free",ncol=ncol)+ylab("Quantity")+xlab("Year")+
  geom_line(aes(x=year,y=obs),linetype="dashed",size=0.7, na.rm=TRUE)
  } else {
    d.$S50 = (d.$sel)
    p=ggplot(d.,aes(x=year,y=data,group=sel))+geom_line(aes(colour=sel), na.rm=TRUE)+
      geom_line(aes(x=year,y=obs),linetype="dashed",size=0.7, na.rm=TRUE)+
      facet_wrap(~qname, scales="free",ncol=ncol)+ylab("Quantity")+xlab("Year")+
      scale_color_discrete(limits = unique(d.$sel))
    }
  
  return(p)
  } #}}}

#{{{
#' mleg()
#
#' function to adjust ggplot legend for multiplots
#'
#' @param size legend.key.size (cm)
#' @param height legend.key.height (cm)
#' @param width legend.key.width (cm)
#' @param titel legend title size
#' @param test legend text size
#' @return legformat   
#' @export

mleg = function(size=0.5,height=0.5,width=0.5,title=6,text=6){
  
  legformat = theme(legend.key.size = unit(size, 'cm'), #change legend key size
           legend.key.height = unit(height, 'cm'), #change legend key height
           legend.key.width = unit(width, 'cm'), #change legend key width
           legend.title = element_text(size=title), #change legend title font size
           legend.text = element_text(size=text)) #change legend text font size

  return(legformat)
}  
#}}}

