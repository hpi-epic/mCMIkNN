##### Experiments on a small synthetic network

rm(list = ls())
source("source.R")

# network definition
generateMixedNetwork = function(n){
  m = 4
  tCov = 0.7
  covM = matrix(c(1,tCov,tCov,1), nrow=2, ncol=2)
  A = rexp(n=n, rate=1)
  B = sample(0:(m-1),n,replace=T)
  C = rep(0,n)
  D = rep(0,n)
  E = rep(0,n)
  F = rep(0,n)
  G = rep(0,n)
  for(i in 1:n){
    C[i] = rbinom(1, prob=0.5, size=B[i])
    D[i] = rnorm(1, mean=B[i], sd=1) - (m / 2.0)
    E[i] = rexp(1, rate=1/(C[i]+1))
    F[i] = D[i]^(round(C[i] / 2)) + rnorm(1)
    dummy = (sign(E[i]-1) + 1) / 2
    if(dummy == 0){
      G[i] = rnorm(1, mean=A[i], sd=1)
    }else{
      G[i] = rpois(1, lambda=A[i]) * dummy
    }
  }
  df = data.frame(A,B,C,D,E,F,G)
  suffStat = list(dm=as.matrix(df), type=c(1,0,0,1,1,1,1), n=n)
  return(suffStat)
}

####
# true edges in the network
from = c("A", "B", "B", "C", "C", "D", "E")
to = c("G", "C", "D", "E", "F", "F", "G")
true_arcs = cbind(from, to)

FF = "results/network_Chisq99.tab"
testI = CMIp.Chisq99 # select test that should be evaluated (see algorithm/CMI_pvals.R)

alpha = 0.01 ## fixed for proposed approach

## init
source("experiments/causal_graph_learning/evaluate_networks.R")

V = c("A", "B", "C", "D", "E", "F", "G")
set.seed(1)
sizes = c(100, 500, 1000, 2000, 5000, 10000)
prec_mean = rep(0, length(sizes))
recall_mean = rep(0, length(sizes))
f1_mean = rep(0, length(sizes))
for(l in 1:length(sizes)){
    numnodes = 7
    print(sizes[l])
    loops = 20
    pLoops = rep(0, loops)
    rLoops = rep(0, loops)
    fLoops = rep(0, loops)
    for(i in 1:loops){
        print(i)
        dat = generateMixedNetwork(sizes[l])
        pc.fit <- pc(suffStat = dat, indepTest = testI, alpha=alpha, labels = V, verbose = TRUE)
        arcs = getEdges(res=pc.fit,nnodes=numnodes,V=V)
        stats = compareToTrueNet(arcs=arcs, trueArcs=true_arcs)
        if(stats$found == 0){
            pLoops[i] = 0
        }else{
            pLoops[i] = stats$tp / stats$found
        }
        rLoops[i] = stats$tp / stats$correct
        fLoops[i] = F1(pLoops[i],rLoops[i])
    }
    f1_mean[l] = mean(fLoops)
    recall_mean[l] = mean(rLoops)
    prec_mean[l] = mean(pLoops)
    df = data.frame(Size=sizes, Prec=prec_mean, Rec=recall_mean, F1=f1_mean)
    write.table(df, file=FF, row.names=F, quote=F, col.names=T, sep="\t")
}



