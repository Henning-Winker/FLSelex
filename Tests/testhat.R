
data(ple4)
ple4_mtf <- stf(ple4, nyears = 10)
# Now the stock goes up to 2018
summary(ple4_mtf)



ple4_sr <- fmle(as.FLSR(ple4, model="bevholt"), control=list(trace=0))
plot(ple4_sr)

f_status_quo <- mean(fbar(ple4)[,as.character(2015:2017)])
f_status_quo

ctrl_target <- data.frame(year = 2018:2027,
                          quant = "f",
                          value = 0.1)

ctrl_f <- fwdControl(ctrl_target)
ctrl_f

a = ple4_f_sq <- fwd(ple4_mtf, control = ctrl_f, sr = ple4_sr)

plot(a)

library(FLCore)
library(FLBRP)
library(ggplot2)
library(FLasher)
library(FLSelex)

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
bt = selex.backtest(pars,stk,Fref="Fmsy")
plotprjselex(bt)



# Deterministic foward projection
fw = selex.fwd(pars,stk,sr=bh)
plotprjselex(fw)


# add ssr function to brp.selex
# SRR
sr = as.FLSR(stk,model=bevholt)
bh = fmle(sr)
plot(FLSRs(bh=bh))

brps.sr = brp.selex(pars,stock=stk,sr=bh)
plotFselex(brps.sr)

ploteqselex(brps.sr)
ploteqselex(brps.sr,panels=4)
