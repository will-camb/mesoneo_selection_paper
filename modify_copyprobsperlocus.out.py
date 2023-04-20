import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobsperlocus_location",
                    help="location of individual copyprobsperlocus.out output from painting",
                    required=True)
args = parser.parse_args()
copyprobsDF = pd.read_csv(args.copyprobsperlocus_location, delim_whitespace=True, low_memory=False)
Hap_start_sites = copyprobsDF.loc[copyprobsDF.pos == 'HAP'].index.tolist()
copyprobsDF.drop('pos', axis=1, inplace=True)
bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
names = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
d = dict(enumerate(names, 1))
mask = ~copyprobsDF.index.isin(Hap_start_sites)
copyprobsDF_valid = copyprobsDF[mask]
cols = copyprobsDF.columns.tolist()
for n, i in enumerate(cols):
    copyprobsDF.loc[mask, cols[n]] = np.vectorize(d.get)(np.digitize(copyprobsDF_valid[cols[n]].astype('float64'), bins))
copyprobsDF.iloc[Hap_start_sites[0]] = str(copyprobsDF.iloc[Hap_start_sites[0]][1]) + "," + str(int(copyprobsDF.iloc[Hap_start_sites[0]][0])) + ","
copyprobsDF.iloc[Hap_start_sites[1]] = str(copyprobsDF.iloc[Hap_start_sites[1]][1]) + "," + str(int(copyprobsDF.iloc[Hap_start_sites[1]][0])) + ","
transposed_copyprobsDF = copyprobsDF.iloc[0:Hap_start_sites[1]].T.append(copyprobsDF.iloc[Hap_start_sites[1]:].reset_index(drop=True).T)
transposed_copyprobsDF.to_csv(args.copyprobsperlocus_location, index=False, header=False, sep=" ")