# ICSL

## Setup and Example

To setup the code execute the following command.
```
git clone git@github.com:JohannesHuegle/ICSL.git
cd src
pip3 install -r requirements.txt
```
You can also install the independence tests via
```
pip3 install .
```
which lets you later call
```
import indeptests
```

Note, within the requirements, we added manm-cs only to allow rerunning some of the code for the experiments and results shown in the associated publication. This is not required for execution of the conditional independence test mCMIkNN.

After setup, you can execute the following command as an example.
`python3 main.py -i ../data/coolingData_1000.csv`

The main.py script takes the following parameter as input:
Parameter | Default | Description
--- | --- | ---
-i | - | Input File location, with column description in header
-l | - | Maximum size of separation set
-c | 2 | Number of cores to execute tests in parallel
-p | 100 | Number of permutations $M_{perm}$
-t | - | Optional transformation of data. Options are "normalize", "rank", "standardize", "uniform" default is none.
-k | 25 | Value $k_{CMI}$ of knn in cmi estimation
-r | 5 | Value $K_{perm}$ of knn in local permutation scheme
-a | 0.01 | Nominal value $\alpha$, i.e., significance level

## Contributor
* [Daniel Thevessen](https://github.com/danthe96)
* [Johannes Huegle](https://github.com/JohannesHuegle)
* [Christopher Hagedorn](https://github.com/ChristopherSchmidt89)
