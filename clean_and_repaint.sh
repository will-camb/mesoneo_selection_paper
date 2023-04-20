#!/bin/bash
cd $PBS_O_WORKDIR
module load anaconda3/4.4.0

chrlist=`seq 1 22`
for chr in $chrlist; do
    # length=`head -n 1 $chr.master_all_copyprobsperlocus.txt`
    length=`awk -F ',' 'NR==2 { print $0 }' phasefiles/$chr.merged.phase`
    awk -F ',' -v nsnps="$length" 'length($3) == nsnps' $chr.master_all_copyprobsperlocus.txt > temp.$chr.master_all_copyprobsperlocus.txt
    cut -f1 -d"," temp.$chr.master_all_copyprobsperlocus.txt | uniq -c > temp.$chr.counts
done

python3 << END
import pandas as pd
import numpy as np

ordered_all_pop_ids_mapped = pd.read_csv("ordered_all_pop_ids_mapped", header=None, sep=" ")
ordered_all_pop_ids_mapped_list = ordered_all_pop_ids_mapped[0].tolist()
to_paint = list()
for i in range(1,23):
    painted = pd.read_csv("temp."+str(i)+".counts", sep=" ", header=None, error_bad_lines=False)[[5,6]]
    painted_successfully_list=painted[painted[5].astype(str)=='20'][6].tolist()
    to_paint_temp = np.setdiff1d(ordered_all_pop_ids_mapped_list, painted_successfully_list)
    to_paint.extend([s for s in to_paint_temp if "UKBB" in s]) # NB only for UKB painting, not ancients!
    print("Done assessment of chr " + str(i) + " and results written to to_paint; now processing chr " + str(i+1))
to_paint = list(set(to_paint))
chunkcounts_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.chunkcounts.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
chunklengths_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.chunklengths.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
regionchunkcounts_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
regionsquaredchunkcounts_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
mutationprobs_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.mutationprobs.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
to_paint.extend([s for s in np.setdiff1d(ordered_all_pop_ids_mapped_list, chunkcounts_list) if "UKBB" in s])
to_paint.extend([s for s in np.setdiff1d(ordered_all_pop_ids_mapped_list, chunklengths_list) if "UKBB" in s])
to_paint.extend([s for s in np.setdiff1d(ordered_all_pop_ids_mapped_list, mutationprobs_list) if "UKBB" in s])
to_paint.extend([s for s in np.setdiff1d(ordered_all_pop_ids_mapped_list, regionchunkcounts_list) if "UKBB" in s])
to_paint.extend([s for s in np.setdiff1d(ordered_all_pop_ids_mapped_list, regionsquaredchunkcounts_list) if "UKBB" in s])
to_paint = list(set(to_paint))
commands = pd.read_csv("paintvspanel1by1_commands.txt", header=None, sep=" ")
commands = commands[commands[2].isin(to_paint)]
commands.to_csv("repaint.paintvspanel1by1_commands.txt", sep=" ", header=False, index=False)
print("New painting commands written to repaint.paintvspanel1by1_commands.txt, now deleting lines in results files containing the individuals to be repainted")

## Now delete lines in results files containing the individuals to be repainted
for i in range(1,23):
  dummy=pd.read_csv("temp." + str(i) + ".master_all_copyprobsperlocus.txt", header=None, error_bad_lines=False)
  dummy = dummy[~dummy[0].isin(to_paint)]
  dummy.to_csv(str(i) + ".master_all_copyprobsperlocus.txt", header=False, index=False)
print("Done for master_all_copyprobsperlocus.txt files, now looking at chunkcounts etc")

chunkcounts = pd.read_csv("ordered_all_pop_ids_mapped.allchr.chunkcounts.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
chunkcounts = chunkcounts[~chunkcounts[0].isin(to_paint)]
chunkcounts.to_csv("ordered_all_pop_ids_mapped.allchr.chunkcounts.out", sep=" ", header=False, index=False)

chunklengths = pd.read_csv("ordered_all_pop_ids_mapped.allchr.chunklengths.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
chunklengths = chunklengths[~chunklengths[0].isin(to_paint)]
chunklengths.to_csv("ordered_all_pop_ids_mapped.allchr.chunklengths.out", sep=" ", header=False, index=False)

regionchunkcounts = pd.read_csv("ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
regionchunkcounts = regionchunkcounts[~regionchunkcounts[0].isin(to_paint)]
regionchunkcounts.to_csv("ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out", sep=" ", header=False, index=False)

regionsquaredchunkcounts = pd.read_csv("ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
regionsquaredchunkcounts = regionsquaredchunkcounts[~regionsquaredchunkcounts[0].isin(to_paint)]
regionsquaredchunkcounts.to_csv("ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out", sep=" ", header=False, index=False)

mutationprobs = pd.read_csv("ordered_all_pop_ids_mapped.allchr.mutationprobs.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
mutationprobs = mutationprobs[~mutationprobs[0].isin(to_paint)]
mutationprobs.to_csv("ordered_all_pop_ids_mapped.allchr.mutationprobs.out", sep=" ", header=False, index=False)

## And for the per-chromosome results
for i in range(1,23):
  chunkcounts = pd.read_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".chunkcounts.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
  chunkcounts = chunkcounts[~chunkcounts[0].isin(to_paint)]
  chunkcounts.to_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".chunkcounts.out", sep=" ", header=False, index=False)
  chunklengths = pd.read_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".chunklengths.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
  chunklengths = chunklengths[~chunklengths[0].isin(to_paint)]
  chunklengths.to_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".chunklengths.out", sep=" ", header=False, index=False)
  regionchunkcounts = pd.read_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".regionchunkcounts.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
  regionchunkcounts = regionchunkcounts[~regionchunkcounts[0].isin(to_paint)]
  regionchunkcounts.to_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".regionchunkcounts.out", sep=" ", header=False, index=False)
  regionsquaredchunkcounts = pd.read_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".regionsquaredchunkcounts.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
  regionsquaredchunkcounts = regionsquaredchunkcounts[~regionsquaredchunkcounts[0].isin(to_paint)]
  regionsquaredchunkcounts.to_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".regionsquaredchunkcounts.out", sep=" ", header=False, index=False)
  mutationprobs = pd.read_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".mutationprobs.out", header=None, sep=" ", error_bad_lines=False, dtype = str)
  mutationprobs = mutationprobs[~mutationprobs[0].isin(to_paint)]
  mutationprobs.to_csv("perchrom_results/ordered_all_pop_ids_mapped.chr"+str(i)+".mutationprobs.out", sep=" ", header=False, index=False)
END

rm -r temp.*

echo "Success! Now run repaint.paintvspanel1by1_commands.txt [e.g. cat repaint.paintvspanel1by1_commands.txt | parallel] or split if large"
