
data("ple4")
stk=ple4
plotselage(stk)
Sa = selage(stk)

fit = fitselex(Sa)
plotselex(fit)
pars = varselex(fit$par,stk,step=0.1,type="shift")
plotselex(pars,stk)
brps = brp.selex(pars,stk)
plotFselex(brps)
ploteqselex(brps)
# SRR
sr = as.FLSR(stk,model=bevholt)
bh = fmle(sr)
plot(FLSRs(bh=bh))

# add ssr function to brp.selex
brps.sr = brp.selex(pars,stock=stk,sr=bh)
ploteqselex(brps.sr)
ploteqselex(brps.sr,panels=4)
