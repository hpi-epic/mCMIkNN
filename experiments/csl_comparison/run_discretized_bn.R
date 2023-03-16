library(bnlearn, quietly = T)
library(dplyr, quietly = T)
library('Ckmeans.1d.dp')
require("parallel")
library('igraph')

### setup directory for learned cgms
dir.create(file.path('./', './results_discretized/'), showWarnings = FALSE)
setwd(file.path('./', './results_discretized/'))

### get all filenames from sample dir
files = list.files('../data_generation/csl_data_normalized/')
for (file in files) {
	df <- read.csv(paste0("../data_generation/csl_data_normalized/",file), header=TRUE, check.names=FALSE, sep=",")
	generate.category.data <- function(X){
        X.categories <- c()
        for (i in c(1:ncol(X))){
            y = X[,i]
            y.categories <- Ckmeans.1d.dp(x=y, k=c(2:min(10,length(unique(y)))), y=1, method="quadratic", estimate.k="BIC")
            X.categories <- cbind(X.categories,y.categories$cluster-1)
        }
        colnames(X.categories) <- colnames(X)
        return(as.data.frame(X.categories))
    }
    cat_df <- generate.category.data(df)
    matrix_df <- cat_df%>%dplyr::mutate_all(funs(if(length(unique(.))<11) as.factor(.)  else as.numeric(as.numeric(.))))
    result = pc.stable(matrix_df, debug=FALSE, test="x2", alpha=0.05, max.sx=Inf)
    write_graph(as.igraph(result), paste0('./results_discretized/',substr(file,1,nchar(file)-4),'.gml'), format = 'gml')
}