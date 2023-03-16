## Experiments

This folder contains all the relevant jupyter notebooks and scripts to recreate the experiment of the paper associated with ```mCMIkNN```.
- Everything concerning data generation is found within the folder ```data_generation```. Note, we only provide example data within this repository, to keep memory footprint low.
- Everything concerning the calibration experiment (Sec. 5.2) is found within the folder ```calibration```.
- The robustness experiment scripts (Sec. 5.2) are found in the folder ```robustness```.
- All scripts to compare the CI tests (Sec. 5.3) are found within the folder ```ci_test_comparison```. Note, you might have to install certain packages in R or python to run the other concidered CI test implementations.
- All scripts relevant for the compiration of the CI tests' performance in the context of causal discovery (Sec. 5.4) are found in the folder ```csl_comparison```. Note, you might have to install certain packages in R or python to run the other concidered CI test implementations in the context of constraint-based causal structure learning. We do not provide any installation scripts.
