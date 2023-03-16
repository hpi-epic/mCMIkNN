## continuous functions
random_fun = function(x,simple=F){
  maxD = 4
  if(simple){
    maxD = 2
  }
  dice = sample.int(4,1)
  if(dice == 1){
    return(x)
  } else if(dice == 2){
    return(x^2)
  } else if(dice == 3){
    return(x^3)
  } else if(dice == 4){
    return(tanh(x))
  } else{
    stop("error: random_fun!")
  }
}

random_constant = function(maxC=4){
  dice1 = runif(1, min=1, max=maxC)
  dice2 = sample.int(2,1) # dice1 or 1/dice1
  dice3 = sample.int(2,1) # * 1 or * (-1)
  const = dice1
  if(dice3 == 2){
    const = const * (-1)
  }
  if(dice2 == 2){
    const = 1 / const
  }
  return(const)
}
simple_constant = function(maxC=4){
  dice1 = runif(1, min=1, max=maxC)
  return(dice1)
}

simple_function = function(x,sd=0.1){
  n = length(x)
  z = random_constant(maxC=2) * random_fun(x) + rnorm(n,sd=sd)
  return(z)
}

random_additive_fun = function(x,y,sd=0.1){
  n = length(x)
  val = random_constant() * random_fun(x) + random_constant() * random_fun(y)
  val = val + rnorm(n,sd=sd)
  return(val)
}

multiplicative_sign_fun = function(x,y){
  n = length(x)
  val = sign(x*y) * rexp(n=n, rate=1/sqrt(2))
  return(val)
}

## mixed functions
mixed_cont_to_bin_sign = function(x,y,noise_level=0.1){
  n = length(x)
  z = sign(x*y)
  noise_ind = sample(1:n, floor(n * noise_level))
  z[noise_ind] = z[noise_ind] * (-1)
  return(z)
}

check_t = function(t,minn){
  if(length(t) < 2){
    return(FALSE)
  }else{
    for(i in 1:length(t)){
      if(t[i] < minn){
        return(FALSE)
      }
    }
    return(TRUE)
  }
}

mixed_cont_to_cat_modulo_simple = function(x,y,noise_level=0.1){
  n = length(x)
  z = (round(simple_constant(maxC=2) * x)) + (1 +  (round(simple_constant(maxC=2) * y))) * 10
  ## randomly assign noise_level percent of the values a random value in the domain
  noise_ind = sample(1:n, floor(n * noise_level))
  domainZ = unique(sort(z))
  noise_val = sample(domainZ, length(noise_ind), replace=T)
  z[noise_ind] = noise_val
  return(z)
}

### its actually categorical not binary
mixed_cont_to_cat_modulo = function(x,y,noise_level=0.1){
  n = length(x)
  z = abs(round(2 * random_constant() * x + random_constant() * y)) %% 4
  t = table(z)
  minn = ceiling(n * 0.1)
  while(!check_t(t=t,minn=minn)){
    z = abs(round(2 * random_constant() * x + random_constant() * y)) %% 4
    t = table(z)
  }
  ## randomly assign noise_level percent of the values a random value in the domain
  noise_ind = sample(1:n, floor(n * noise_level))
  domainZ = unique(sort(z))
  noise_val = sample(domainZ, length(noise_ind), replace=T)
  z[noise_ind] = noise_val
  return(z)
}

xor_scaling = function(x,y,sd=0.1){
  n = length(x)
  zdash = (x + y) %% 2
  z = (rpois(n=1,lambda=5) + 1) * zdash * rnorm(n, sd=sd)
  
}

xor_shift = function(x,y,sd=0.1){
  n = length(x)
  zdash = (x + y) %% 2
  z = (rpois(n=1,lambda=5) + 1) * zdash + rnorm(n, sd=sd)
  
}

mean_shift_fun = function(x,sd=0.1){
  n = length(x)
  val = random_constant() + random_fun(x) + rnorm(n,sd=sd)
  return(val)
}

generate_normal_or_uniform = function(n){
  dice = sample.int(2,1)
  if(dice == 1){
    return(rnorm(n))
  }else{
    return(runif(n,min=-2,max=2))
  }
}

# Fig.2 arxiv
test1 = function(n, rate=0.1){
  # X -> Z -> Y  ----> I(X,Y | Z) = 0
  # x = exp(10) --- 0.5 rate is ok and much faster
  # z = pois(x)
  # y = binom(z,0.5)
  x = rexp(n, rate=rate)
  z = rep(0,n)
  y = rep(0,n)
  for(i in 1:n){
    z[i] = rpois(1, lambda=x[i])
    y[i] = rbinom(1, p=0.5, size=z[i])
  }
  return(data.frame(x,y,z))
}

# fig.3 arxiv
test2 = function(n,m=3){
  x = sample(0:(m-1),n,replace=T)
  y = rep(0,n)
  z = rbinom(n, p=0.5, size=3)
  for(i in 1:n){
    y[i] = runif(1, min=x[i], max=x[i]+2)
  }
  return(data.frame(x,y,z))
}
## zero inflated poisson (neurips paper)
test3 = function(n,p=0.15){
  x = rexp(n, rate=1)
  y = rep(0,n)
  z = rbinom(n, p=0.5, size=3)
  for(i in 1:n){
    coin = runif(1, min=0, max=1)
    if(coin < p){
      y[i] = 0
    }else{
      y[i] = rpois(1, lambda=x[i])
    }
  }
  return(data.frame(x,y,z))
}

# fig.5 arxiv
test5 = function(n, corV=0.8, mixtureProb = c(1/2,1/2), discrete_prob=c(0.8,0.2)){
  z = rbinom(n, size = 3, prob = 0.2)
  mixture_index = sample(c(1,2), size = n, replace = T, prob = mixtureProb)
  
  xy = matrix(rep(0,n*2),nrow=2)
  
  xy[,mixture_index==1] = mvrnorm(n = sum(mixture_index==1), mu = c(0,0),
                                  Sigma = matrix(c(1,corV,corV,1),ncol=2))
  
  xy[,mixture_index==2] = expand.grid(c(1,-1),
                                      c(1,-1))[sample(
                                        1:4,size=sum(mixture_index==2),replace = T,prob=c(0.4,0.1,0.1,0.4)),] %>% t
  x = xy[1,]
  y = xy[2,]
  return(data.frame(x,y,z))
}
