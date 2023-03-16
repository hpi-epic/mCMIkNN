## Causal Discovery - Comparison (Sec. 5.4)

This folder contains the script to learn causal graphical models (CGM) using constraint-based causal discovery and applying different CI tests, that have been used in the CI test comparison and are relevant for the evaluation in our paper.
The learned CGMs obtained with a distinct CI test are stored separately. The evaluation script computes error measures comparing the ground truth undirected skeleton graph with the learned undirected skeleton graph. Based upon the error measures the evaluation script compiles the figures provided within the paper and its appendix.

Note that we use a slightly adpated version for KCIT, as we encountered an issue for some extreme cases.
In detail we encountered the following issue:
```Error in eivy %*% diag(sqrt(eig_Kyz)) : non-conformable arguments```
Explained:
Error occured as ```eivy``` turned to be a single column and is no longer stored as a matrix, as a result ```sqrt(eig_Kyz)``` is no vector, but a single value and ```diag()``` returns a matrix determined by the size of ```sqrt()``` with diag set to ```1```. Hence, we wrote custom script that copies the functions of KCIT and included a check for both ```eivy``` and ```eivx```, i.e., ```if (is.matrix())``` then run normal, ```else``` do not perform matrix multiplication but normal multiplication and transform result via ```as.matrix()```. We also added try catch around the execution to obtain further issues, due to our data.
Hence, other erroneous executions will be discarded in the evaluation.


