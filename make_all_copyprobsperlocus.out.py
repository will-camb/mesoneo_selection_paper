import pandas as pd 
import argparse 
import os

parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="which chromosome",
                    required=True)
parser.add_argument("-phasefile",
                    help="phasefile to read number of SNPs",
                    required=True)
parser.add_argument("-o",
                    help="location to save all_copyprobsperlocus.out to",
                    required=True)
args = parser.parse_args()
file_name = args.chr + '.all_copyprobsperlocus.txt'
if not os.path.isfile(os.path.join(args.o, file_name)):
    all_copyprobsperlocusDF = pd.DataFrame(columns=pd.read_csv(args.phasefile, nrows=2, header=None).iloc[1])
    all_copyprobsperlocusDF.to_csv(os.path.join(args.o, file_name), sep=' ', index=False)
