rm(list = ls())
source("source.R")
n = 1000
x = rnorm(n,0,1)
y = rnorm(n,0,1)
z = sign(x + y + rnorm(n,0,0.1))
data = cbind(x,y,z)

# 1D adaptive histogram
oned_hist = oned_hist_iterative(x,eps=0.05,Kmax = 20)
hist(x,breaks = oned_hist[,1:2] %>% unlist() %>% unique())


# 2D adaptive histogram, iterative
twod_hist_iter = iterative_cmi_greedy_flexible_parallel(cbind(x,y),cores = 2)
plot(x,y)
abline(v=twod_hist_iter[,1:2]%>%unlist%>%unique,
                         h=twod_hist_iter[,3:4]%>%unlist%>%unique,col="red")

# 3D adaptive histogram, iterative
hist_iter_xyz = iterative_cmi_greedy_flexible_parallel(data,isCat = 3,cores=2)

plot(x,y)

plot(x[z==1],y[z==1]);abline(v=hist_iter_xyz[,1:2]%>%unlist%>%unique,
                                     h=hist_iter_xyz[,3:4]%>%unlist%>%unique,col="red")
plot(x[z==-1],y[z==-1]);abline(v=hist_iter_xyz[,1:2]%>%unlist%>%unique,
                                       h=hist_iter_xyz[,3:4]%>%unlist%>%unique,col="red")

#### Estimate CMI for mixture data

source("../experiments/syntheticGenerator.R")
## generate mixture data (see paper Experiment V)
data = test5(n)

## calculate CMI based on multidimensional histograms
## Estimate contains a list of values, where 'lh' is the empirical estimate based on the calculated histogram and the remaining values are different correction terms which we explain in the paper cited in the main README file.
estimates = CMI.estimates(data=data, xind=1,yind=2,zinds=c(3),isCat=c(3), logE=T)

### true CMI estimate for the given data
trueCMI = 0.4 * log(0.4/(0.5^2)) + 0.1 * log(0.1/(0.5^2)) - 0.25 * log((1-(0.8^2)))
### empirical estimate based on histograms
estimates$LH / n
### empitical estimate with chi-squared correction with alpha=0.05
(estimates$LH + estimates$Chisq95) / n
### or
estimate.Chisq95 = CMI.Chisq95(data=data, xind=1,yind=2,zinds=c(3),isCat=c(3), logE=T) / n

## For different correction criteria see file 'algorithm/CMI_estimates.R'

### To use one of the corrected CMI estimates for the PC algorithm from the pcalg package, instead use CMIp.Chisq95 or similar (see 'algorithm/CMI_pvals.R'). This generates a pseudo p-value > 0.01 iff our estimate is <= 0, i.e., indicates independence and a p-value <= 0.01 otherwise. For an example see "experiments/causal_graph_learning/syntheticNonLinear.R".


