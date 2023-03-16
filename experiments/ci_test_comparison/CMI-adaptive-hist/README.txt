The code provided in this package supports the research paper below.

@inproceedings{marx:21:myl,
    title={{Estimating Conditional Mutual Information for Discrete-Continuous Mixtures using Multidimensional Adaptive Histograms}},
    author={Marx, Alexander and Yang, Lincen and van Leeuwen, Matthijs},
    year={2021},
    booktitle = {SDM, Virtual Conference},
    publisher = {SIAM},
}

The code for the algorithm is in ./algorithm. The file ./algorithm/main.R contains some simple examples for learning adaptive histogram from data and estimate CMI from mixture data.
- CMI_estimates.R contains the functions, which compute a CMI estimate based on adaptive histogram models using different correction criteria
- CMI_pvals.R contains wrapper functions, which compute pseudo p-values form estimated CMI_values

The experiments folder contains:
- evaluations of CMI estimates on ground-truth data ("CMI_estimation", Figure 1 & Figure 2);
- causal graph structure learning ("causal_graph_learning", Figure 4)
- evaluations of our estimator as conditional independence test ("conditionalIndep_testing", Figure 6);
To run the experiments, stay in this folder SDM2021/code/ with your console and run, e.g.,
Rscript experiments/CMI_estimation/test_distrib.R
The results will be stored into the 'results' folder. Currently, each test related to independence testing uses the chi-squared correction. This can be changed manually in the corresponding folder.

It is required to have R version >= 3.6.0. You may also need to install additional R packages (see requirements in 'source.R'). If some of the packages are not available through "install.packages", try installing via bioconductor: e.g. https://www.bioconductor.org/packages/release/bioc/html/graph.html



