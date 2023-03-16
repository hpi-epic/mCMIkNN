library(RCIT)
require("parallel")

header <- data.frame("cgmid","samples","sepsetsize","discretenoderatio","method","pvalue","hasedge")
write.table(header, file = './KCIT.csv', append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")

files = list.files('../data_generation/ci_data_normalized/')
cores <- 2

comp <- function(file) {
	elements <- c()
	splitted = strsplit(file,"_")
	cgmid = strtoi(splitted[[1]][1])
	dnr = as.numeric(paste(splitted[[1]][3],'.',splitted[[1]][4],sep=""))
	s = strtoi(splitted[[1]][6])
	hasedge = if (substr(splitted[[1]][7],1,nchar(splitted[[1]][7])) == 'withEdge') 'True' else 'False'
	sample = strtoi(substr(splitted[[1]][8],1,nchar(splitted[[1]][8])-4))
	df <- read.csv(paste0("../data_generation/ci_data_normalized/",file), header=TRUE, check.names=FALSE, sep=",")
		pval <- -1
		tryCatch({
			pval = KCIT(as.vector(t(df[paste(s)])),as.vector(t(df[paste(s+1)])),as.matrix(df[paste(seq(0,s-1))]))
		}, error=function(e){
		### cannot use CI test, keep pval as -1, handle it in evaluation script as an error case
		})
		element <- data.frame(cgmid,sample,s,dnr,'KCIT',pval,hasedge)
		elements <- rbind(elements, element)

	return(elements)
}

results <- mclapply(files, comp, mc.cores=cores)
results_cleaned <- c()
for (e in results) {
	results_cleaned <- rbind(results_cleaned, e)
}
write.table(results_cleaned, file = './KCIT.csv', append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep=",")
