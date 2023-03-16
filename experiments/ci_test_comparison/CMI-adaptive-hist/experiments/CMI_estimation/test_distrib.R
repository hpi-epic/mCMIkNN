#### Experiments I-V in the paper

rm(list = ls())
source("source.R")
source("experiments/syntheticGenerator.R")
steps = 1:10
loops = 1:100
fact = 100

## Test 1: uniform dependence
mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
m = 5
trueI = log(m) - (m-1)*log(2)/m
index = 1
set.seed(1)
## end modified
for(s in steps){
    print(s)
    n = s * fact
    estimates = data.frame(LH=rep(0, length(loops)), JIC=rep(0, length(loops)), qNML=rep(0, length(loops)), fNML=rep(0, length(loops)), Chisq99=rep(0, length(loops)), Chisq95=rep(0, length(loops)))
    for(l in loops){
        dd = test2(n,m=m)
        estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(1),logE=T)
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
        for(i in 2:6){
          estimates[l,i] = max((estimate[[1]] + estimate[[i]]) / n,0)
        }
    }
    for(i in 1:6){
      mean_vec[index,i+1] = mean(estimates[,i])
      sd_vec[index,i+1] = sd(estimates[,i])
      meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / length(loops)
    }
    print(mean_vec[index,])
    print(sd_vec[index,])
    print(meanse_vec[index,])
    index = index + 1
}
write.table(mean_vec, file="results/test1_mean5.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/test1_sd5.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/test1_meanse5.tab", sep="\t", row.names=F, col.names=T, quote=F)

### test 2 independent
### independent scenario.. chain
mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))

trueI = 0
index = 1
set.seed(1)
## end modified
for(s in steps){
    print(s)
    n = s * fact
    estimates = data.frame(LH=rep(0, length(loops)), JIC=rep(0, length(loops)), qNML=rep(0, length(loops)), fNML=rep(0, length(loops)), Chisq99=rep(0, length(loops)), Chisq95=rep(0, length(loops)))
    for(l in loops){
      dd = test1(n,rate=0.5)
        estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(2,3), logE=T)
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
        for(i in 2:6){
          estimates[l,i] = max((estimate[[1]] + estimate[[i]]) / n, 0)
        }
    }
    for(i in 1:6){
      mean_vec[index,i+1] = mean(estimates[,i])
      sd_vec[index,i+1] = sd(estimates[,i])
      meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / length(loops)
    }
    print(mean_vec[index,])
    print(sd_vec[index,])
    print(meanse_vec[index,])
    index = index + 1
}
write.table(mean_vec, file="results/test2_mean.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/test2_sd.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/test2_meanse.tab", sep="\t", row.names=F, col.names=T, quote=F)

### test 3 independent
### zero inflated poisson
mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))

p = 0.15
trueI = (1-p) * 0.3012
index = 1
set.seed(1)
## end modified
for(s in steps){
    print(s)
    n = s * fact
    estimates = data.frame(LH=rep(0, length(loops)), JIC=rep(0, length(loops)), qNML=rep(0, length(loops)), fNML=rep(0, length(loops)), Chisq99=rep(0, length(loops)), Chisq95=rep(0, length(loops)))
    for(l in loops){
      dd = test3(n,p=p)
        estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(2), logE=T)
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
        for(i in 2:6){
          estimates[l,i] = max((estimate[[1]] + estimate[[i]]) / n,0)
        }
    }
    for(i in 1:6){
      mean_vec[index,i+1] = mean(estimates[,i])
      sd_vec[index,i+1] = sd(estimates[,i])
      meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / length(loops)
    }
    print(mean_vec[index,])
    print(sd_vec[index,])
    print(meanse_vec[index,])
    index = index + 1
}
write.table(mean_vec, file="results/test3_mean15.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/test3_sd15.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/test3_meanse15.tab", sep="\t", row.names=F, col.names=T, quote=F)


### test 4 two dim Gaussian with covariance r
mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))

tCov = 0.6
trueI = -0.5 * log(1-tCov^2)
covM = matrix(c(1,tCov,tCov,1), nrow=2, ncol=2)
index = 1
require("MASS")
set.seed(1)
## end modified
for(s in steps){
    print(s)
    n = s * fact
    estimates = data.frame(LH=rep(0, length(loops)), JIC=rep(0, length(loops)), qNML=rep(0, length(loops)), fNML=rep(0, length(loops)), Chisq99=rep(0, length(loops)), Chisq95=rep(0, length(loops)))
    for(l in loops){
        dd = mvrnorm(n=n, Sigma=covM, mu=c(0,0))
        estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(), logE=T)
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
        for(i in 2:6){
          estimates[l,i] = max((estimate[[1]] + estimate[[i]]) / n,0)
        }
    }
    for(i in 1:6){
      mean_vec[index,i+1] = mean(estimates[,i])
      sd_vec[index,i+1] = sd(estimates[,i])
      meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / length(loops)
    }
    print(mean_vec[index,])
    print(sd_vec[index,])
    print(meanse_vec[index,])
    index = index + 1
}
write.table(mean_vec, file="results/test4_mean6.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/test4_sd6.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/test4_meanse6.tab", sep="\t", row.names=F, col.names=T, quote=F)

### test 5 independent
### mixture
mean_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
sd_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))
meanse_vec = data.frame(Samples=steps*fact, LH=rep(0, length(steps)), JIC=rep(0, length(steps)), qNML=rep(0, length(steps)), fNML=rep(0, length(steps)), Chisq99=rep(0, length(steps)), Chisq95=rep(0, length(steps)))

trueI = 0.4 * log(0.4/(0.5^2)) + 0.1 * log(0.1/(0.5^2)) - 0.25 * log((1-(0.8^2)))
index = 1
set.seed(1)
## end modified
for(s in steps){
    print(s)
    n = s * fact
    estimates = data.frame(LH=rep(0, length(loops)), JIC=rep(0, length(loops)), qNML=rep(0, length(loops)), fNML=rep(0, length(loops)), Chisq99=rep(0, length(loops)), Chisq95=rep(0, length(loops)))
    for(l in loops){
      dd = test5(n)
        estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(3),isCat=c(3), logE=T)
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
        for(i in 2:6){
          estimates[l,i] = max((estimate[[1]] + estimate[[i]]) / n,0)
        }
    }
    for(i in 1:6){
      mean_vec[index,i+1] = mean(estimates[,i])
      sd_vec[index,i+1] = sd(estimates[,i])
      meanse_vec[index,i+1] = sum((estimates[,i] - trueI)^2) / length(loops)
    }
    print(mean_vec[index,])
    print(sd_vec[index,])
    print(meanse_vec[index,])
    index = index + 1
}
write.table(mean_vec, file="results/test5_mean15.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/test5_sd15.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/test5_meanse15.tab", sep="\t", row.names=F, col.names=T, quote=F)
