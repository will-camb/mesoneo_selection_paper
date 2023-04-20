#  Script designed for Neale Lab GWAS sumstats
#  This script can LD prune, and looks at variants in all individuals (not just homozygous as in _v4)


import pandas as pd
import argparse
import os


def analyse_anc(merged_phase_copyprobs_temp, anc, chrom, pval, iteration, phenotype_file):
    global list_of_SNPs
    skipped_snps = 0
    phase_sum = merged_phase_copyprobs_temp.loc[:, (slice(None), 'phase')].sum().tolist()
    phase_sum_dict = {}
    for n, i in enumerate(list_of_SNPs):
        phase_sum_dict[i] = phase_sum[n]
    alt_mapping = {}
    output = pd.DataFrame(index=list_of_SNPs, columns=['alt', 'ref'])
    for i in list_of_SNPs:
        try:
            minor_allele = GWAS_variants.loc[GWAS_variants['pos'] == i].minor_allele.item()
        except ValueError:  # This skips SNPs with more then one entry in the GWAS file
            print("ValueError for minor_allele for SNP at position " + str(i) + " for " + str(anc) + " chr" +
                  str(chrom) + "! Skipping this SNP (probably because >1 entry in the GWAS file)")
            skipped_snps += 1
            continue
        alt_allele = GWAS_variants.loc[GWAS_variants['pos'] == i].alt.item()
        #  Check if the minor allele = alt allele
        if minor_allele == alt_allele:
            #  minor is alt
            #  Find if minor is 1 or 0 in phase
            if phase_sum_dict[i] > (len(copyprobs_haps) / 2):
                alt_mapping[i] = 0
            else:
                #  major is alt
                alt_mapping[i] = 1
        elif minor_allele != alt_allele:
            #  major is alt
            if phase_sum_dict[i] > (len(copyprobs_haps) / 2):
                alt_mapping[i] = 1
            else:
                alt_mapping[i] = 0
        alt_sum = merged_phase_copyprobs_temp.loc[merged_phase_copyprobs_temp[i]['phase'] == alt_mapping[i]].loc[:,
                  (i, 'copyprobs')].sum()
        ref_sum = merged_phase_copyprobs_temp.loc[merged_phase_copyprobs_temp[i]['phase'] != alt_mapping[i]].loc[:,
                  (i, 'copyprobs')].sum()
        output.at[i, 'alt'] = alt_sum
        output.at[i, 'ref'] = ref_sum
    output['maf'] = output['alt'] / (output['alt'] + output['ref'])
    output['maf'] = output['maf'].fillna(0)
    output = output.reset_index().rename(columns={'index': 'pos'})
    maf_GWAS = pd.merge(output, GWAS_variants, on='pos')
    maf_GWAS['maf_x_beta'] = maf_GWAS['maf'] * maf_GWAS['beta']
    PRS = maf_GWAS['maf_x_beta'].sum()
    number_of_SNPs = len(list_of_SNPs)
    results_list.append(
        [args.copyprobs_file, phenotype_file, anc, chrom, PRS, number_of_SNPs, skipped_snps, pval, iteration,
         maf_GWAS['maf'].tolist(), maf_GWAS['beta'].tolist()])


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file",
                    help="Should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile",
                    help="Should be named in form chr#.merged.phase.gz",
                    required=True)
parser.add_argument("-idfile",
                    help="Directory of idfile used in chromopainter ('individual pop 1')",
                    required=True)
parser.add_argument("-phenotypes",
                    help="File with list of phenotypes being looked at; NB extension must be .gz and not.bgz!",
                    required=True)
parser.add_argument("-chr",
                    help="Chromosome number",
                    required=True)
parser.add_argument("-anc",
                    help="Ancestry being analysed",
                    required=True)
parser.add_argument("-bootstrap",
                    help="Whether to bootstrap or not; can take two values: True or False",
                    required=True)
args = parser.parse_args()
if args.bootstrap == 'True':
    bootstrap = True
else:
    bootstrap = False

# Read in copyprobs
col_names = pd.read_csv(str(args.copyprobs_file), sep=" ", nrows=0).columns
types_dict = {'0': str}
types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
anc_copyprobs = pd.read_csv(str(args.copyprobs_file), sep=" ", dtype=types_dict)
anc_copyprobs.set_index("0", inplace=True)
anc_copyprobs.columns = anc_copyprobs.columns.tolist()[::-1]  # Reverse column names because they are mis-labelled;
#  so when correct, column names are descending
copyprobs_haps = list()
for h in range(int(anc_copyprobs.shape[0] / 2)):
    copyprobs_haps.extend([1, 2])
print("*** Successfully loaded copyprobs file for " + args.anc + " chr" + str(args.chr) + "***")
# Read in idfile
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
# Read in phasefile
phase = pd.read_csv(str(args.phasefile), header=None, sep=" ", dtype='int8')
phase.columns = anc_copyprobs.columns.tolist()[::-1]  # Phasefile cols should be ascending
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
phase_haps = list()
for h in range(int(phase.shape[0] / 2)):
    phase_haps.extend([1, 2])
# Load variants file to find which is alt allele (not always minor!)
variants = pd.read_csv("/willerslev/datasets/UKBiobank/NealeV2/variants.tsv.gz",
                       sep='\t',
                       usecols=['variant', 'chr', 'alt'],
                       dtype={'variant': 'string', 'chr': 'string', 'alt': 'string'})
variants = variants.loc[variants['chr'] == str(args.chr)]
ldetect = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/PRS_calculation/ldetect-data/EUR/fourier_ls-chr"
                      + str(args.chr) + ".bed", sep='\t')
ldetect.rename(columns=lambda x: x.strip(), inplace=True)
results_list = list()
# LOOP through GWAS files
phenotypes = pd.read_csv(args.phenotypes, header=None, dtype=str)[0].tolist()
phenotypes = [s.strip() for s in phenotypes]
for file in phenotypes:
    print("***Looking at " + str(file) + "***")
    GWAS = pd.read_csv(file,
                       sep='\t',
                       usecols=['variant', 'beta', 'pval', 'minor_allele'],
                       dtype={'variant': 'string', 'beta': 'float', 'pval': 'float', 'minor_allele': 'string'})
    GWAS[['chr', 'pos', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
    GWAS = GWAS.loc[GWAS['chr'] == str(args.chr)]
    GWAS = GWAS.loc[GWAS['pos'].isin(phase.columns.tolist())]  # NB about 80% of painted SNPs are in GWAS file
    GWAS_variants = pd.merge(variants[['variant', 'alt']], GWAS, on='variant')
    # for p in [0.0001, 0.001, 0.01]:
    for p in [5e-8]:
        GWAS_variants_pval = GWAS_variants.loc[GWAS_variants['pval'] < p]
        if GWAS_variants_pval.empty:
            print("No SNPs pass p-val threshold for " + str(args.chr))
            continue
        best_per_block = pd.DataFrame(columns=GWAS_variants.columns)
        for index, row in ldetect.iterrows():
            block = GWAS_variants_pval.loc[(GWAS_variants_pval['pos'].astype(int) >= row['start']) & (
                    GWAS_variants_pval['pos'].astype(int) < row['stop'])].reset_index()
            try:
                block = block.reset_index()
                # best_per_block.loc[index] = block.iloc[block['beta'].abs().idxmax()]  # Take highest beta
                best_per_block.loc[index] = block.iloc[block['pval'].idxmin()]  # Take lowest p-val
            except ValueError:  # For when there are no SNPs that pass p-val threshold in the block
                continue
        list_of_SNPs = best_per_block['pos'].drop_duplicates().tolist()
        phase_temp = phase[list_of_SNPs]
        phase_temp = phase_temp.reset_index().rename(columns={'index': 'ID'})
        phase_temp['haps'] = phase_haps
        anc_copyprobs_temp = anc_copyprobs[list_of_SNPs]
        anc_copyprobs_temp = anc_copyprobs_temp.reset_index().rename(columns={'0': 'ID'})
        anc_copyprobs_temp['haps'] = copyprobs_haps
        merged_phase_copyprobs = pd.merge(phase_temp, anc_copyprobs_temp, on=['ID', 'haps'])
        reorder = [[str(x) + "_x", str(x) + "_y"] for x in list_of_SNPs]
        reorder = [item for sublist in reorder for item in sublist]
        merged_phase_copyprobs = merged_phase_copyprobs[['ID', 'haps'] + reorder]
        iterables = [["ID"] + list_of_SNPs, ["phase", "copyprobs"]]
        merged_phase_copyprobs.columns = pd.MultiIndex.from_product(iterables, names=['first', 'second'])
        merged_phase_copyprobs.set_index('ID', inplace=True)
        if bootstrap:
            for bs in range(1000):
                temp = merged_phase_copyprobs.sample(n=merged_phase_copyprobs.shape[0], replace=True)
                analyse_anc(temp, str(args.anc), args.chr, p, bs, file)
        else:
            analyse_anc(merged_phase_copyprobs, str(args.anc), args.chr, p, bootstrap, file)

print("***Success! Now writing results to output file***")
if not os.path.exists("PRS_calculations_bootstrapped_v3_" + str(args.phenotypes)):
    open("PRS_calculations_bootstrapped_v3_" + str(args.phenotypes), 'a').close()
PRS_calculations = pd.DataFrame.from_records(results_list)
PRS_calculations.to_csv("PRS_calculations_bootstrapped_v3_" + str(args.phenotypes), mode='a', header=False, index=False,
                        sep=" ")
