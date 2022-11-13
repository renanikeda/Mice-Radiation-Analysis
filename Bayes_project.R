library("coda")
library("rjags")
dat = read.table("http://users.stat.ufl.edu/~winner/data/micerad.dat", header=F)
dat$V2 <- dat$V2+1
dat$V4 <- (dat$V1-mean(dat$V1))/sd(dat$V1)

head(dat)
table(dat$V3,dat$V1)
table(dat$V3,dat$V2)
lmod = lm(dat$V3~dat$V1+dat$V2)
anova(lmod)
plot(jitter(resid(lmod)))
y_hat = predict(lmod,dat[,1:2])
decision0 = as.numeric(y_hat > 0.5)

plot(jitter(dat$V3)~jitter(dat$V4),xlab="neutrons",ylab="died")
plot(jitter(dat$V3)~jitter(dat$V2),xlab="therapy",ylab="died")

table(dat$V3,dat$V2)

mod1_string = " model {
  for (i in 1:length(died)) {
    died[i] ~ dbern(p[i])
    logit(p[i]) = b0 + b[1]*neutrons[i] + b[2]*treatment[i]
  } 
  
  b0 ~ dnorm(0.0, 1.0/1e4)
  
  for (j in 1:2) {
    b[j] ~ ddexp(0.0, sqrt(2.0)) 
  }
  
} "
data1_jags = list(neutrons = dat$V4, treatment = dat$V2, died = dat$V3)

params1 = c("b0","b")

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                       variable.names=params1,
                       n.iter=1e4)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))
windows()
plot(mod1_sim)
gelman.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)
(dic1 = dic.samples(mod1, n.iter=1e3))

pm_params1 = colMeans(mod1_csim)
logdied1_hat = pm_params1["b0"] + pm_params1["b[1]"]*dat$V4 +pm_params1["b[2]"]*dat$V2
died1_hat = 1/(1+exp(-logdied1_hat))
decision1 = as.numeric(died1_hat>0.5)
resid1 = data1_jags$died - decision1
plot(jitter(resid1))
SQ1 = sum((data1_jags$died - decision1)^2)


mod2_string = " model {
  for (i in 1:length(died)) {
    died[i] ~ dbern(p[i])
    logit(p[i]) = int[treatment[i]] + b[1]*neutrons[i] 
  } 
  
  for (k in 1:max(treatment)) {
    int[k] ~ dnorm(mu, prec)
  }
  
  b[1] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
  
  mu ~ dnorm(0, 1.0/1e4)
  prec ~ dgamma(1.0,1.0)
  sig = 1/sqrt(prec)

} "
data2_jags = list(neutrons = dat$V4, treatment = dat$V2, died = dat$V3)

params2 = c("int","b","mu","sig")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=1e4)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))
windows()
plot(mod2_sim,ask=T)
gelman.diag(mod2_sim)
effectiveSize(mod2_sim)
autocorr.plot(mod2_sim)
(dic2 = dic.samples(mod2, n.iter=1e3))

(mean(mod2_csim[,"int[1]"] > mod2_csim[,"int[2]"]))

pm_params2 = colMeans(mod2_csim)
logdied2_hat = pm_params2[2:3][dat$V2] + pm_params2["b"]*dat$V4
died2_hat = 1/(1+exp(-logdied2_hat))
decision2 = as.numeric(died2_hat>0.5)
resid2 = data2_jags$died - died2_hat
plot(jitter(resid2))
SQ2 = sum((data2_jags$died - decision2)^2)

test = c(mean(dat$V4),1,1)
logtest = mod2_csim[,2]+ mod2_csim[,"b"]*mean(dat$V4)
result = 1/(1+exp(-logtest))
mean(result>0.5)