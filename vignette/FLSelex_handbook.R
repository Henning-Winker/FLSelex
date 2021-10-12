#-----------------------------------------------------
# FLSelex: Test Script FLSelex Handbook
# Henning Winker JCR-EC
# Email: henning.winker@ec.europa.eu
# Licence: e EUPL 1.1.
#-----------------------------------------------------



## ---- eval=FALSE------------------------------------------------------------------------------------------
## installed.packages("devtools")
## 
## devtools::install_github("flr/FLCore")
## 
## devtools::install_github("flr/FLBRP")
## 
## devtools::install_github("flr/FLasher")
## 
## devtools::install_github("flr/ggplotFL")
## 
## devtools::install_github("henning-winker/FLSelex")
## 


## ---------------------------------------------------------------------------------------------------------
library(FLCore)
library(FLBRP)
library(FLasher)
library(FLSelex)
library(ggplotFL)


## ----fig1, fig.height=3.5, fig.cap = "Observed $S_a$ = $F_a$/$F_{max}$ for North Sea place over the recent 5 years"----

data(ple4)

plotselage(ple4,nyears=5)



## ---------------------------------------------------------------------------------------------------------
Sa = selage(ple4,nyears=3)
Sa



## ----fig2, fig.height=3.5, fig.cap = "Observed $S_a$ and fitted $S_a$ for North Sea place using `fitselex()`"----

Sa = selage(ple4,nyears=3)

fit = fitselex(Sa)

plotselex(sel=fit,Sa=Sa)



## ----fig3, fig.height=3.5, fig.cap = "Observed $S_a$ and fitted $S_a$ for North Sea plaice illustrating the compounds of the piece-wise `selex` function. logis: logistic of acending curve ($S_{50}$, $S_{95}$), hnorm: unadjusted halfnormal ($S_{max}$, $D_{cv}$) for the descending curve, height:  adjusted height of the halfnormal ($D_{min}$)"----

plotselex(sel=fit,Sa=Sa,compounds=TRUE)



## ----fig4, fig.height=3.5, fig.cap = "Observed $S_a$ and fitted $S_a$ for North Sea place, assuming a simplified logistic selectivity"----

ogivefit = as.ogive(fit)

plotselex(sel=ogivefit,Sa=Sa)



## ----fig5, fig.height=3.5, fig.cap = "Cranking the ascending slope of the estimated  selectivity curve by varying $S_{50}$ "----

crank = varselex(pars=fit$par,stock=ple4,step=0.1,type="crank")

plotselex(sel=crank,Sa=Sa)



## ----fig6, fig.height=3.5, fig.cap = "Shifting the estimated  selectivity curve in its unchanged shape  by varying $S_{50}$, $S_{95}$ and $S_{max}$"----

shift = varselex(pars=fit$par,stock=ple4,step=0.1,type="shift")

plotselex(sel=shift,Sa=Sa)



## ----fig7, fig.height=3.5, fig.cap = "Shifting the estimated  selectivity curve in its unchanged shape  by varying $S_{50}$, $S_{95}$ and $S_{max}$"----

dyn = varselex(pars=fit$par,stock=ple4,step=0.1,type="dynamic")

plotselex(sel=dyn,Sa=Sa)



## ---------------------------------------------------------------------------------------------------------

range(ple4)[c("minfbar","maxfbar")]
test = fbar2f(ple4)
range(test)[c("minfbar","maxfbar")]



## ---------------------------------------------------------------------------------------------------------

pars = varselex(fit$par,ple4,type="dynamic")

brps = brp.selex(pars,ple4)

class(brps)


## ----fig8, fig.height=5.5, fig.cap = "Plots showing the trade-offs between $S_50$ and $F$ with respect to relative yield-per-recruit and the  spawning ratio potential (SPR, or spawning biomass per recruit). Dashed lines connecting at solid black dots denote the expected outcome of current $F$ and $S_a$ at equilibrium"----

ploteqselex(brps)



## ----fig9, fig.height=3.5, fig.cap = "Fit of a Beverholt Model to the spawing stock biomass and recruitment estimates for North Sea plaice",message=FALSE,warning=FALSE,result='hide'----

sr = as.FLSR(ple4,model=bevholt)
bh = fmle(sr)
plot(FLSRs(bh))+theme(legend.position = "right")



## ----fig10, fig.height=4.5, fig.cap = "Plots showing the trade-offs between $S_50$ and $F$ with respect to relative yield and and $SSB$ based on a Beverton-Holt SSR. Dashed lines connecting at solid black dots denote the expected outcome of current $F$ and $S_a$ at equilibrium"----

brps.sr = brp.selex(pars,ple4,sr=bh)
ploteqselex(brps.sr)



## ----fig12, fig.height=6.5,fig.width=5, fig.cap = "Plots showing relative changes in (top) $YPR$ and $SPR$ and (bottom) yield and $SSB$ over range of $S50$ values under the current $F$ "----

p1=plotFselex(brps,what="Fref")+ggtitle("Per-Recruit")
p2=plotFselex(brps.sr,what="Fref")+ggtitle("Beverton-Holt SRR")
gridExtra::grid.arrange(p1,p2,ncol=1)



## ----fig14, fig.height=5.5, fig.cap = "Plots showing the responses of Catch, Harvest Rate, Percentage of juveniles in the catch and SSB to changes in selectivity for determistic future forecasts over 30 years"----

bt = selex.fwd(pars,ple4,sr=bh,fyears=30)
plotprjselex(bt)


## ----fig16, fig.height=5.5, fig.cap = "Plots showing the responses of Catch, Harvest Rate, Percentage of juveniles in the catch and SSB to changes in selectivity if they had actually been implemeted in 2009"----

bt = selex.backtest(pars,ple4,byears=10)
plotprjselex(bt)

