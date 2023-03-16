#### Experiment VI in the paper

rm(list = ls())
source("source.R")
source("experiments/syntheticGenerator.R")
loops = 1:100
n = 2000
nz = 0:4

## Test 1: uniform dependence
mean_vec = data.frame(NZ=nz, LH=rep(0, length(nz)), JIC=rep(0, length(nz)), qNML=rep(0, length(nz)), fNML=rep(0, length(nz)), Chisq99=rep(0, length(nz)), Chisq95=rep(0, length(nz)))
sd_vec = data.frame(NZ=nz, LH=rep(0, length(nz)), JIC=rep(0, length(nz)), qNML=rep(0, length(nz)), fNML=rep(0, length(nz)), Chisq99=rep(0, length(nz)), Chisq95=rep(0, length(nz)))
meanse_vec = data.frame(NZ=nz, LH=rep(0, length(nz)), JIC=rep(0, length(nz)), qNML=rep(0, length(nz)), fNML=rep(0, length(nz)), Chisq99=rep(0, length(nz)), Chisq95=rep(0, length(nz)))
m = 5
trueI = log(m) - (m-1)*log(2)/m
index = 1
set.seed(1)
## end modified
for(s in nz){
    print(s)
    estimates = data.frame(LH=rep(0, length(loops)), JIC=rep(0, length(loops)), qNML=rep(0, length(loops)), fNML=rep(0, length(loops)), Chisq99=rep(0, length(loops)), Chisq95=rep(0, length(loops)))
    for(l in loops){
        dd = test2(n,m=m)
        if(s > 1){
          for(gz in 2:s){
            dd = data.frame(dd, z=rbinom(n, p=0.5, size=3))
          }
        }
        estimate = NULL
        if(s == 0){
          estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(),isCat=c(1),logE=T)
        }else{
          estimate = CMI.estimates(data=dd, xind=1,yind=2,zinds=c(1:s + 2),isCat=c(1,1:s + 2),logE=T)
        }
        estimates[l,1] = estimate[[1]] / n ## pure likelihood
        for(i in 2:6){
          estimates[l,i] = (estimate[[1]] + estimate[[i]]) / n
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
write.table(mean_vec, file="results/testmd_mean5_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(sd_vec, file="results/testmd_sd5_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
write.table(meanse_vec, file="results/testmd_meanse5_n2000.tab", sep="\t", row.names=F, col.names=T, quote=F)
