# This function load all the packages and functions needed for the discretization. 
require("SCCI")
require("parallel")
require(dplyr)

Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")


source("generate_candidate_cuts.R")
source("refine_candidate_cuts_exclude_empty_bins.R")
source("modified_log.R")
source("iterative_cmi_greedy_flexible_parallel.R")

source("multi_hist_splitting_seed_based_simple.R")
source("oned_hist_iterative.R")
source("CMI_estimates.R")
source("CMI_pvals.R")
source("utils.R")
