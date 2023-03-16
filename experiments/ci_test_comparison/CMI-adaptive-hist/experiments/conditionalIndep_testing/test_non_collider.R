#### Test different generating mechanisms that for a non-collider structure

rm(list = ls())
source("source.R")
source("experiments/syntheticGenerator.R")
# modify this:
fileF="results/test_non_collider_chisq99.tab"
testI = CMIp.Chisq99 # select test that should be evaluated (see algorithm/CMI_pval.R)

alpha = 0.01 ## fixed for proposed approach
steps = 1:10
loops = 1:100
fact = 100

## Test 1: random continuous chain
acc = rep(0, length(steps))
index = 1
set.seed(1)
## could be modified
noise = 0.1  ## standard deviation of additive noise
## end modified
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = generate_normal_or_uniform(n)
        y = random_constant() * random_fun(x) + rnorm(n,sd=noise)
        z = random_constant() * random_fun(y) + rnorm(n,sd=noise)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(1,1,1), n=n)
        pval = testI(x=1,y=3,S=c(2),suffStat=dat)
        if(pval >= alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
### results df
df = data.frame(Samples=steps*fact, ContChain=acc, ContCC=rep(0, length(steps)), MixedChain=rep(0, length(steps)), MixedCC=rep(0, length(steps)))
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)


### Test difficult collider, marginal independencies...
acc = rep(0, length(steps))
index = 1
set.seed(1)
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = generate_normal_or_uniform(n)
        y = random_constant() * random_fun(x) + rnorm(n,sd=noise)
        z = random_constant() * random_fun(x) + rnorm(n,sd=noise)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(1,1,1), n=n)
        pval = testI(x=2,y=3,S=c(1),suffStat=dat)
        if(pval >= alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
df$ContCC = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

### X, Y continuous to binary via sign
acc = rep(0, length(steps))
index = 1
set.seed(1)
## could be modified
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        dd = test1(n,rate=0.5)
        dat = list(dm=dd, type=c(1,0,0), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval >= alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
df$MixedChain = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

### X,Y continuous to Z categorical via modulo
acc = rep(0, length(steps))
index = 1
m = 5
set.seed(1)
## could be modified
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = sample(0:(m-1),n,replace=T)
        y = rep(0,n)
        z = rep(0,n)
        const = random_constant()
        for(i in 1:n){
          y[i] = runif(1, min=x[i], max=x[i]+2)
          z[i] = rnorm(1, sd=x[i]) + const
        }
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(0,1,1), n=n)
        pval = testI(x=2,y=3,S=c(1),suffStat=dat)
        if(pval >= alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
df$MixedCC = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)
