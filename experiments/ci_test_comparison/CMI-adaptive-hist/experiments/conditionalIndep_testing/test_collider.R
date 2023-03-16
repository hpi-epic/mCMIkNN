#### Test different generating mechanisms that form a collider structure

rm(list = ls())
source("source.R")
source("experiments/syntheticGenerator.R")
# modify this:
fileF="results/test_collider_chisq99.tab"
testI = CMIp.Chisq99 # select test that should be evaluated (see algorithm/CMI_pval.R)

alpha = 0.01 ## fixed for proposed approach
# the test is assumed to have a wrapper such as CMIp.XX, that is get a suffStat
#end modify
steps = 1:10
loops = 1:100
fact = 100

## Test 1: non-linear continuous relationships
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
        y = generate_normal_or_uniform(n)
        z = random_additive_fun(x,y,sd=noise)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(1,1,1), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval < alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
### results df
df = data.frame(Samples=steps*fact, Poly=acc, Sign=rep(0, length(steps)), SignM=rep(0, length(steps)), ModM=rep(0, length(steps)), XorM=rep(0, length(steps)), XorSd=rep(0, length(steps)))
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
        x = rnorm(n)
        y = rnorm(n)
        z = multiplicative_sign_fun(x,y)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(1,1,1), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval < alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
df$Sign = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

### X, Y continuous to binary via sign
acc = rep(0, length(steps))
index = 1
set.seed(1)
## could be modified
noise = 0.1   ### percentage of randomly assigned values
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = rnorm(n)
        y = rnorm(n)
        z = mixed_cont_to_bin_sign(x,y,noise_level=noise)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(1,1,0), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval < alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
df$SignM = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

### X,Y continuous to Z categorical via modulo
acc = rep(0, length(steps))
index = 1
set.seed(1)
## could be modified
noise = 0.1   ### percentage of randomly assigned values
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = rnorm(n)
        y = rpois(n, lambda=sample(2,1)) + 1
        z = round(x %% y)
        noise_ind = sample(1:n, floor(n * noise))
        domainZ = unique(sort(z))
        noise_val = sample(domainZ, length(noise_ind), replace=T)
        z[noise_ind] = noise_val
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(1,0,0), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval < alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}
df$ModM = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

### X,Y binary, xor determines mean shift
acc = rep(0, length(steps))
index = 1
set.seed(1)
## could be modified
noise = 1  ## standard deviation of z that is shifted
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = round(runif(n, min=0, max=1))
        y = round(runif(n, min=0, max=1))
        z = xor_shift(x,y,sd=noise)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(0,0,1), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval < alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}

df$XorM = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

### X,Y binary, xor determines sd scaling
acc = rep(0, length(steps))
index = 1
set.seed(1)
## could be modified
noise = 0.1  ## standard deviation of z that is scaled
for(s in steps){
    print(s)
    n = s * fact
    corr_count = 0
    for(l in loops){
        x = round(runif(n, min=0, max=1))
        y = round(runif(n, min=0, max=1))
        z = xor_scaling(x,y,sd=noise)
        dd = as.matrix(data.frame(x,y,z))
        dat = list(dm=dd, type=c(0,0,1), n=n)
        pval = testI(x=1,y=2,S=c(3),suffStat=dat)
        if(pval < alpha){
            corr_count = corr_count + 1
        }
    }
    acc[index] = corr_count / length(loops)
    print(acc[index])
    index = index + 1
}

df$XorSd = acc
write.table(df, file=fileF, sep="\t", row.names=F, col.names=T, quote=F)

