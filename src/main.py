import argparse
import pandas as pd

from csl.pc_parallel import parallel_stable_pc
from mCMIkNN import mCMIkNN

parser = argparse.ArgumentParser("Parallel PC Algorithm on Mixed Discrete-Continuous Data")
parser.add_argument('-a', '--alpha', help='Signficance Level used for CI Test', default=0.01, type=float)
parser.add_argument('-l', '--level', help='Maximum Level of the Run', default=None, type=int)
parser.add_argument('-c', '--process_count', help='Number of Processes used in the Run - Main process excluded', default=2, type=int)
parser.add_argument('-p', '--permutations', help='Number of Permutations used for a CI Test', default=100, type=int)
parser.add_argument('-t', '--transform', help='Choice of transform are "normalize", "rank", "standardize", "uniform" default is none', default=None, choices=['normalize', 'rank', 'standardize', 'uniform'], type=str)
parser.add_argument('-k', '--kcmi', help='KNN during cmi estimation.', type=int, default=25)
parser.add_argument('-r', '--kperm', help='KNN during restricted local permutation computation.', type=int, default=5)
required = parser.add_argument_group('required arguments')
required.add_argument('-i', '--input_file', help='Input File Name', required=True)


if __name__ == '__main__':

	args = parser.parse_args()

	input_data = pd.read_csv(args.input_file)
	citest = mCMIkNN(args.kcmi, args.kperm, args.permutations, args.transform)

	res = parallel_stable_pc(input_data, citest, args.alpha, args.process_count, args.level)
