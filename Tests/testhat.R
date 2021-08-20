
library(FLCore)
library(FLBRP)
library(ggplot2)
library(FLasher)

data("ple4")
stk=ple4
plotselage(stk)
Sa = selage(stk)

fit = fitselex(Sa)
plotselex(fit)
pars = varselex(fit$par,stk,step=0.2,type="dynamic")
plotselex(pars,stk)
brps = brp.selex(pars,stk)
plotFselex(brps)
ploteqselex(brps)
# backtest
bt = selex.backtest(pars,stk,sr=bh,Fref="Fmsy")
plotprjselex(bt)
# Deterministic foward projection
fw = selex.fwd(pars,stk)
plotprjselex(fw)


# SRR
sr = as.FLSR(stk,model=bevholt)
bh = fmle(sr)
plot(FLSRs(bh=bh))

# add ssr function to brp.selex
brps.sr = brp.selex(pars,stock=stk,sr=bh)
plotFselex(brps.sr)

ploteqselex(brps.sr)
ploteqselex(brps.sr,panels=4)
