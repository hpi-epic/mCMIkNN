# todo fix manual import of cmipchisq95
library(R.utils)
require("parallel")
require("SCCI")
require(dplyr)
library("MASS")
Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")
source("./CMI-adaptive-hist/algorithm/generate_candidate_cuts.R")
source("./CMI-adaptive-hist/algorithm/refine_candidate_cuts_exclude_empty_bins.R")
source("./CMI-adaptive-hist/algorithm/modified_log.R")
source("./CMI-adaptive-hist/algorithm/iterative_cmi_greedy_flexible_parallel.R")

source("./CMI-adaptive-hist/algorithm/multi_hist_splitting_seed_based_simple.R")
source("./CMI-adaptive-hist/algorithm/oned_hist_iterative.R")
source("./CMI-adaptive-hist/algorithm/CMI_estimates.R")
source("./CMI-adaptive-hist/algorithm/CMI_pvals.R")
source("./CMI-adaptive-hist/algorithm/utils.R")


header <- data.frame("cgmid","samples","sepsetsize","discretenoderatio","method","pvalue","hasedge")
write.table(header, file = './cmipchisq95.csv', append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")

files = list.files(path='.../data_generation/ci_data_normalized/')
cores <- 2

comp <- function(file) {
	elements <- c()
	splitted = strsplit(file,"_")
	cgmid = strtoi(splitted[[1]][1])
	dnr = as.numeric(paste('0.',splitted[[1]][4],sep=""))
	s = strtoi(splitted[[1]][6])
	hasedge = if (substr(splitted[[1]][7],1,nchar(splitted[[1]][7])) == 'withEdge') 'True' else 'False'
	sample = strtoi(substr(splitted[[1]][8],1,nchar(splitted[[1]][8])-4))
	df <- read.csv(paste0("../data_generation/ci_data_normalized/",file), header=TRUE, check.names=FALSE, sep=",")
		matrix_df <- data.matrix(df)
	    type <- rep(1, ncol(matrix_df))
	    for(i in 1:ncol(matrix_df)) {
	        if(length(unique(matrix_df[ , i])) < 11) {
	            type[i] = 0
	        }
	    }
	    sufficient_stats = list(dm=matrix_df, type=type, n=nrow(df))
		pval <- -1
		tryCatch({
			pval <- withTimeout( {CMIp.Chisq95(s+1,s+2,seq(1,s),sufficient_stats)}, timeout = 600)
		}, error=function(e){
			### cannot use CI test, keep pval as -1, handle it in evaluation script as an error case
			pval <- -1
		})
		element <- data.frame(cgmid,sample,s,dnr,'cmipchisq95',pval,hasedge)
		elements <- rbind(elements, element)
	return(elements)
}

results <- mclapply(files, comp, mc.cores=cores)
results_cleaned <- c()
for (e in results) {
	results_cleaned <- rbind(results_cleaned, e)
}
write.table(results_cleaned, file = './cmipchisq95.csv', append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")
