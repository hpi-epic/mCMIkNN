library(micd, quietly = T)
require("parallel")
library(pcalg, quietly = T)
library(dplyr, quietly = T)
library(igraph)

### setup directory for learned cgms
dir.create(file.path('./', './results_MICD/'), showWarnings = FALSE)
setwd(file.path('./', './results_MICD/'))

### get all filenames from sample dir
files = list.files('../data_generation/csl_data_normalized/')
for (file in files) {
	df <- read.csv(paste0("../data_generation/csl_data_normalized/",file), header=TRUE, check.names=FALSE, sep=",")
	matrix_df <- df%>%dplyr::mutate_all(funs(if(length(unique(.))<10) as.factor(.)  else as.numeric(.)))
	result = pc(suffStat=matrix_df, verbose=FALSE,
            indepTest=mixCItest, m.max=Inf,
            p=ncol(df), alpha=0.05, numCores=1, skel.method="stable")
    write_graph(graph_from_graphnel(getGraph(result), name = TRUE), paste0('./results_MICD/',substr(file,1,nchar(file)-4),'.gml'), format = 'gml')
}