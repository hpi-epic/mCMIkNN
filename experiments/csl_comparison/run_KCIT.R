require("parallel")
library(pcalg, quietly = T)
library(igraph)


RCITnormalize <- function(mat){

  if (is.null(nrow(mat))){mat = matrix(mat);}

  mat = apply(mat, 2, function(x) if (sd(x)>0){(x - mean(x)) / sd(x)} else{x-mean(x);})


}



U_KCI <-function(x,y){

if (sd(x)>0) x=RCITnormalize(x);
if (sd(y)>0) y=RCITnormalize(y);


T=length(x);

if (length(x)>500){
  width = sqrt(2)*median(as.vector(dist(cbind(x[1:500],y[1:500]))));
} else {
  width = sqrt(2)*median(as.vector(dist(cbind(x,y))));
}

theta = 1/(width^2);

H =  diag(T) - (1/T)*matrix(1,T,T);

Kx = RBF_kernel(x, theta); Kx = H %*% Kx %*% H;
Ky = RBF_kernel(y, theta); Ky = H %*% Ky %*% H;

Sta = sum(diag(Kx %*% Ky));

mean_appr = sum(diag(Kx)) * sum(diag(Ky)) /T;
var_appr = 2* sum(diag(Kx %*% Kx)) %*% sum(diag(Ky %*% Ky))/T^2;
k_appr = (mean_appr^2)/var_appr;
theta_appr = var_appr/mean_appr;


p_val = 1-pgamma(Sta, shape=k_appr, scale=theta_appr);

return(p_val)

}

RBF_kernel <- function(x,width){

  n2=(dist(cbind(x,x),diag = TRUE, upper = TRUE))^2;
  n2=as.matrix(n2);

  width=width/2;
  kx = exp(-n2*width);

  return(kx)

}


KCIT <- function(x,y,z,Bootstrap=TRUE){

  if( length(z)==0 || missing(z) ){
    p_val = U_KCI(x,y);
    #print(p_val)
    return(p_val)
  } else{

    x=RCITnormalize(x);
    y=RCITnormalize(y);
    z=RCITnormalize(z);

    T=length(x);

    lambda = 1E-3;
    Thresh = 1E-5;

    if (T <= 200) {
      #width = 0.8;
      width = 1.2;
    } else if (T > 1200){
      #width = 0.3;
      width = 0.7;
    }else{
      #width = 0.5;
      width = 0.4;
      }

    if (length(x)>500){
    width = median(as.vector(dist(cbind(x[1:500],y[1:500]))));
    } else {
    width = median(as.vector(dist(cbind(x,y))));
    }

    if (is.null(dim(z)[2])) {
    	D = 1;
    } else{
    	D = dim(z)[2];
    }

    theta = 1/(width^2 * D);

    H =  diag(T) - (1/T)*matrix(1,T,T);

    Kx = RBF_kernel(cbind(x, z/2), theta); Kx = H %*% Kx %*% H;

    Ky = RBF_kernel(y, theta); Ky = H %*% Ky %*% H;

    Kz = RBF_kernel(z, theta); Kz = H %*% Kz %*% H;

    #P1 = (diag(T)-Kz %*% chol2inv(chol(Kz + lambda*diag(T))));
    P1 = (diag(T)-Kz %*% solve(Kz + lambda*diag(T)));
    Kxz = P1 %*% Kx %*% t(P1);
    Kyz = P1 %*% Ky %*% t(P1);

    Sta = sum(diag(Kxz %*% Kyz));

    df = sum(diag(diag(T)-P1));

    listxz = eigen((Kxz+t(Kxz))/2,symmetric=TRUE);
    eig_Kxz=listxz[[1]]; eivx=listxz[[2]]

    listyz = eigen((Kyz+t(Kyz))/2,symmetric=TRUE);
    eig_Kyz=listyz[[1]]; eivy=listyz[[2]]


    IIx = which(eig_Kxz > max(eig_Kxz) * Thresh);
    IIy = which(eig_Kyz > max(eig_Kyz) * Thresh);
    eig_Kxz = eig_Kxz[IIx];
    eivx = eivx[,IIx];
    eig_Kyz = eig_Kyz[IIy];
    eivy = eivy[,IIy];

    if (is.matrix(eivx)) {
    	eiv_prodx = eivx %*% diag(sqrt(eig_Kxz));
    } else {
    	eiv_prodx = as.matrix(eivx * sqrt(eig_Kyx));
    }
    if (is.matrix(eivy)) {
    	eiv_prody = eivy %*% diag(sqrt(eig_Kyz));
    } else {
    	eiv_prody = as.matrix(eivy * sqrt(eig_Kyz));
    }

    Num_eigx = dim(eiv_prodx)[2];
    Num_eigy = dim(eiv_prody)[2];
    Size_u = Num_eigx * Num_eigy;

    uu = matrix(0,T,Size_u);


    for (i in 1:Num_eigx){
      for (j in 1:Num_eigy){
        uu[,(i-1)*Num_eigy + j] = eiv_prodx[,i] * eiv_prody[,j];
      }
    }

    if (Size_u > T){
      uu_prod = uu %*% t(uu);
    } else {
      uu_prod = t(uu) %*% uu;
    }

    if (Bootstrap){
      T_BS=5000;
      IF_unbiased=TRUE;

      list_uu = eigen(uu_prod);
      eig_uu =list_uu[[1]];
      II_f = which(eig_uu > max(eig_uu) * Thresh);
      eig_uu = eig_uu[II_f];

      if (length(eig_uu)*T < 1E6){
        f_rand1 = matrix(rnorm(length(eig_uu)*T_BS)^2,length(eig_uu),T_BS);
        Null_dstr = t(eig_uu) %*% f_rand1;

      } else {

        Null_dstr = matrix(0,1,T_BS);
        Length = max(c(floor(1E6/T),100));
        Itmax = floor(length(eig_uu)/Length);
        for (iter in 1:Itmax){
          f_rand1 = matrix(rnorm(Length*T_BS)^2,Length,T_BS);
          Null_dstr = Null_dstr + t(eig_uu[((iter-1)*Length+1):(iter*Length)]) %*% f_rand1;

        }
      }

      sort_Null_dstr = sort(Null_dstr);
      #Cri = sort_Null_dstr[ceiling((1-alpha)*T_BS)];
      p_val = sum(Null_dstr>Sta)/T_BS;
      #print(p_val)
      return(p_val)

    } else {
      mean_appr = sum(diag(uu_prod));
      var_appr = 2*sum(diag(uu_prod^2));
      k_appr = mean_appr^2/var_appr;
      theta_appr = var_appr/mean_appr;
      p_val = 1-pgamma(Sta, shape=k_appr, scale=theta_appr);
      #print(p_val)
      return(p_val)

    }

  }

}


kcitWrapper <- function(x,y,z, suffStat) {
	if (identical(z, integer(0))) {
		U_KCI(as.vector(t(suffStat$dm[x])),as.vector(t(suffStat$dm[y])))
	} else {
		KCIT(suffStat$dm[x],suffStat$dm[y],suffStat$dm[z])
	}
}


### setup directory for learned cgms
dir.create(file.path('./', './results_KCIT/'), showWarnings = FALSE)
setwd(file.path('./', './results_KCIT/'))

### get all filenames from sample dir
files = list.files('../data_generation/csl_data_normalized/')
### set cores
cores <- 2

comp <- function(file) {
  	df <- read.csv(paste0("../data_generation/csl_data_normalized/",file), header=TRUE, check.names=FALSE, sep=",")
  	sufficient_stats <- list(dm=df)
    result = 1
    tryCatch({
    result = pc(suffStat=sufficient_stats, verbose=FALSE,
          indepTest=kcitWrapper, m.max=Inf,
                p=ncol(df), alpha=0.05, numCores=1, skel.method="stable")
    }, error=function(e){
    print(paste0("Error for ", file))
  })
  if (is.numeric(result)) {
    print(paste0(i," ",file))
  } else {
      write_graph(graph_from_graphnel(getGraph(result), name = TRUE), paste0('./results_KCIT/',substr(file,1,nchar(file)-4),'.gml'), format = 'gml')
  }
}

### execute in parallel
x <- mclapply(files, comp, mc.cores=cores)
