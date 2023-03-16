library(bnlearn, quietly = T)
library(dplyr, quietly = T)
library("Ckmeans.1d.dp")

### write header
header <- data.frame("cgmid","samples","sepsetsize","discretenoderatio","method","pvalue","hasedge")
write.table(header, file = './discretized_x2.csv', append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")

files = list.files('../data_generation/ci_data_normalized/')
for (file in files){
	### process each CGM from a unique file
	# get parameters
	splitted = strsplit(file,"_")
	cgmid = strtoi(splitted[[1]][1])
	dnr = as.numeric(paste(splitted[[1]][3],'.',splitted[[1]][4],sep=""))
	s = strtoi(splitted[[1]][6])
	hasedge = if (substr(splitted[[1]][7],1,nchar(splitted[[1]][7])) == 'withEdge') 'True' else 'False'
	sample = strtoi(substr(splitted[[1]][8],1,nchar(splitted[[1]][8])-4))
	df <- read.csv(paste0("../data_generation/ci_data_normalized/",file), header=TRUE, check.names=FALSE, sep=",")
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
		pval <- -1
		tryCatch({
			pval = ci.test(paste(s),paste(s+1),paste(seq(0,s-1)),matrix_df, 'x2')$p.value
		}, error=function(e){
		### cannot use CI test, keep pval as -1, handle it in evaluation script as an error case
		})
		element <- data.frame(cgmid,sample,s,dnr,'discretized-x2',pval,hasedge)
		write.table(element, file = './discretized_x2.csv', append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")
}





