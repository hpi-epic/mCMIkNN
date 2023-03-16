## DATA_GENERATION and DATA

Folder to store the generated data for the experiments, which includes the folders with datasets generated using ```manm-cs``` for experiments concerning the CI tests (Sec. 5.2 and 5.3) and for experiments concerning causal discovery (Sec. 5.4).
The data generated is stored in ```ci_data``` and ```csl_data```, correspondingly. We transform the data using a min-max normalization (see corresponding transformation scripts within this folder) and use the transformed data within our experimental evaluation.
Note, for the data generated for causal discovery, we also store the ground truth graph in the folder ```csl_graph```.

To recreate the experiments, please execute the data generation notebooks provided within this folder. Including the data transformation cells that perform the normalization.

We provide one example data file for each folder.

Structure:

|- ci_data
|- ci_data_normalized
|- csl_data
|- csl_data_normalized
|- csl_graph
|- generate_ci_data.ipynb
|- generate_csl_data.ipynb
