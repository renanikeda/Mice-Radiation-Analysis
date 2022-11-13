install.packages("rjags")
library("rjags")
library("coda")

mod_string = "model{
  for (i in 1:n){
    y[i] ~ dnorm(mu,1.0/sig2) #the 1/sig is the precision, or inv of variance
  }
  mu ~ dt(0,1/1,1)
  sig2 = 1
}"

set.seed(50)
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
n = length(y)

data_jags = list(y = y, n = n)
params = c("mu")
inits = function(){
  inits = list("mu" = 0)
}

mod = jags.model(textConnection(mod_string), data_jags, inits)

update(mod,500)
mod_sim = coda.samples(mod, params, 1e3)

plot(mod_sim)
summary(mod_sim)
