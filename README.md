# mesoneo_selection_paper
Code used in Irving-Pease et al 2023


Painting Manual:

Painting a modern dataset using a pre-defined reference panel

This document outlines the steps necessary to paint a modern dataset using a pre-defined set of reference samples grouped into reference populations. It does not detail how to group the reference samples in the first place. 

Starting point: imputed, phased data in VCF format. It is important that the phasing of the modern data is to a similar standard as the ancient data (and ideally the UKB data too), as poor phasing will mess up the results. 

We can split the process into two main steps:
1.	Data preparation
a.	Merging of datasets
b.	Subsetting and filtering SNPs
c.	File conversions
2.	Painting 
a.	Running painting scripts to record local and global ancestry estimates

Data preparation:

1.	Preparation of modern ‘target’ data:

For the UKB data I used hard genotype calls which were extracted from .bgen files. The main thing to be aware of is to use only software that keeps the phasing of the data (e.g. not certain plink commands…). 

2.	Merge the ancient and modern datasets, subset for selected SNPs

I did this in one step using qctool (https://www.well.ox.ac.uk/~gav/qctool_v2/), merging the ancient and modern data (both in VCF format), and subsetting the positions. For each chromosome, something like:

qctool -g /willerslev/ukbiobank/nonbritish/merged_vcfs/ukbb.haplotype_chrN.GT.chrinfo.samples.vcf.gz -g /willerslev/users-shared/science-snm-willerslev-wl4sn3/step3_postprocessing/step5_release/15062020/chrN.haplotypes.filtered.N1664.D15062020.GLIMPSE.vcf.gz -og merged.N.vcf -os merged.N.samples -compare-variants-by position -flip-to-match-cohort1 -incl-samples keep_samples -incl-positions /willerslev/ukbiobank/SNP_selection/qctool_format/N_SNPs_formatted.txt

where keep_samples lists all ancient and modern IDs that you want to include, and N_SNPs_formatted.txt is a file containing SNPs in the UKB array. This is about ~820k SNPs chosen for markers of specific interest, rare coding variants, and genome-wide coverage. When combined with the imputed ancient data we end up with about 550K SNPs. Chromopainter will not accept missing data.

3.	Filter for MAF > 0.01 and INFO

I filtered for MAF > 0.01 using bcftools. Something like:

bcftools view -q 0.01:minor merged.N.vcf.gz -o merged.N.filtered.vcf.gz

4.	Convert VCF to phase format

Phase format is the chromopainter input format (more info here https://people.maths.bris.ac.uk/~madjl/finestructure/manual.html). There are a few ways to convert a vcf to phasefile format, but they’re all slow (i.e. takes many days/weeks) with large numbers of samples and can’t be parallelised. I found that the best way to do this was to split the vcf by sample, then convert these individual vcf files to phase format and concatenate the results. It’s a little hacky but speeds it up massively. To split the vcf for each chromosome N:

	bcftools plugin split merged.N.filtered.vcf.gz --samples-file samples0 -Oz -o haps.N

Then using script vcf_to_phase.sh (which uses pbwt [https://github.com/richarddurbin/pbwt] and the script impute2chromopainter.pl which is included as part of fineSTRUCTURE), I wrote separate commands to convert each vcf to a phasefile. Using python (pandas) to generate these commands, where file samples contains all sample IDs (modern and ancient):

>>> for chr in range(1,23):
...   df=pd.read_csv("samples", header=None)
...     df_temp=df
...     df_temp[1]="bash vcf_to_phase.sh " + df_temp[0].astype(str) + " " + str(chr)
...     df_temp=df_temp[1]
...     df_temp.to_csv(str(chr)+".commands", header=False, index=False)

Which creates a separate command file for each chromosome which can then be run in parallel. To then concatenate the files in the correct order for chromosome N:

	while read p; do cat phase.N/$p.reduced.phase >> merged.N.phase; done <samples

To get the top three lines for the phase format, we run the same commands as in vcf_to_phase.sh just for one sample but keep the top three lines instead. For example, for chromosome N and sample X:

pbwt -readVcfGT haps.N/X.vcf.gz -writeImputeHapsG haps.N/X.haps
impute2chromopainter.pl haps.N/X.haps temp.N.X.phase
head -n 3 temp.N.X.phase > temp.N.header.phase
cat merged.N.phase >> temp.N.header.phase
mv temp.N.header.phase N.merged.phase

We now have a phasefile with all the samples in the same order as in file samples. The only thing we need to do is edit the number of haplotypes (i.e. 2x number of individuals) which is specified in the first line of the phasefile. Depending on the number, something like:

nhaps=SOME NUMBER
for chr in $chrlist ; do sed -i "1s/.*/$nhaps/" $chr.merged.phase; done

5.	Make recombination input file

We need to make a file specifying recombination rates between our included SNPs. I used the recombination maps from the International Hapmap Consortium phase II release, and the command for each chromosome looked like: 

convertrecfile.pl -M hapmap phasefiles/$chr.merged.phase 2011-01_phaseII_B37/genetic_map_GRCh37_chr$chr.txt recombfiles/$chr.recombfile

convertrecfile.pl is a script included in fineSTRUCTURE.

6.	Make idfile input file

The idfile needs to contain the individuals included in the phasefile in the same order, one per line. The format is:
<NAME> <POPULATION> <INCLUSION> <ignored extra info> 

Where <NAME> and <POPULATION> are strings and <INCLUSION> is 1 to include an individual and 0 to exclude them.  
EXAMPLE IDFILE: 
Ind1 Pop1 1 
Ind2 Pop1 1 
Ind3 Pop2 0 
Ind4 Pop2 1 
Ind5 Pop2 1 

It is best if individuals are named as e.g. UKB1, UKB2, UKB3 etc, otherwise there are some downstream scripts that won’t work. So you need a mapping between the old IDs and the new IDs in this input file. I will send the idfile for the ancient samples but you need to make sure it is in the right order. The first two lines of my idfile looked like:
UKBB1 UKBB 1
UKBB2 UKBB 1

NB the idfile must be in the same order as the phasefile. 

Painting
When the data preparation steps above are done, you should have:
•	A phasefile for each chromosome
•	A recombination file for each chromosome
•	An idfile specifying all your ancient and modern individuals, and which population they are in. 
 
PAINTING PROCESS
Make sure the idfile and phasefile have the target panel at the top and the reference panel below. 

Step 1: Run fineSTRUCTURE on the reference panel only
fs ref_panel.cp -idfile ordered_all_pop_ids_mapped -phasefiles phasefiles/{1..22}.merged.phase -recombfiles ../recombfiles/{1..22}.recombfile -hpc 1 -go

This gives Ne and mu estimates, which can be found in ref_panel.cp file.

Step 2: Run stage03
Edit will_ancientvsancient.cp (including popneinf and popmuinf from fineSTRUCTURE run on ref panel above) and will_ancientvsancient.donors. Also update will_modernvsancient.cp with these numbers. You may want to change the name of the target panel from UKBB to something else, in which case this needs to be reflected in will_modernvsancient.donors and also the idfile. 

Paint ancient panel vs ancient panel to get priors. This paints individuals in reference panel against the panel (but not including itself) using uniform prior first, to obtain the overall copying averages; then repaints using previous run as prior.

Commands:
bash will_03-paintpanelvspanel.sh
This runs will_paint_withinpanel.sh, which in turn runs will_paintsample_withinpanel.sh on each sample.

Step 3: Run stage04
I had a memory problem when running painting in parallel because each thread has to read in the entire phasefiles, so I decided to run in separate batches of 24k. Within each 24k batch, I ran will_3.5-paintvspanel1by1.py which writes commands for 24k individuals, then splits these into separate batchcommands in batch_files folder. To generate commands for this in python:

for j in range(1,410000,24000):
    print("dir=split_"+str(j)+"-"+str((j-1)+24000))
    print("phaselinenumber="+str((j+1)*2))
    print("phaselinenumber2="+str(2*(j+24000)+1))
    print("nhaps=48636")
    print("mkdir $dir")
    print("mkdir $dir/phasefiles")
    print("for chr in $chrlist ; do")
    print("  touch $dir/phasefiles/$chr.merged.phase")
    print("  head -n 3 phasefiles/$chr.merged.phase > $dir/phasefiles/$chr.merged.phase")
    print('  awk "NR>=$phaselinenumber && NR<=$phaselinenumber2" phasefiles/$chr.merged.phase >> $dir/phasefiles/$chr.merged.phase')
    print("  tail -n 636 phasefiles/$chr.merged.phase >> $dir/phasefiles/$chr.merged.phase")
    print('  sed -i "1s/.*/$nhaps/" $dir/phasefiles/$chr.merged.phase')
    print('  echo "Copied lines to $dir/phasefiles/$chr.merged.phase"')
    print("done")
    print('awk "NR>='+str(j)+' && NR<='+str((j-1)+24000)+'" ordered_all_pop_ids_mapped > $dir/ordered_all_pop_ids_mapped')
    print('tail -n 318 ordered_all_pop_ids_mapped >> $dir/ordered_all_pop_ids_mapped')
    print("cp chr1..22_1000inds_test_noflag/cp_panel_scripts/* $dir") 
    print("cp -r recombfiles/ $dir") 
    print("cd $dir")
    print("python3 will_3.5-paintvspanel1by1.py -idfile ordered_all_pop_ids_mapped -cp_panel_scripts ../chr1..22_1000inds_test_noflag/cp_panel_scripts/")
    print("mkdir batch_files")
    print("split -l 1000 -d --additional-suffix=.txt paintvspanel1by1_commands.txt batch_commands")
    print("mv batch_commands* batch_files")
    print("cd ../")

Each of these command files (found in batch_commands, 1000 commands in each) was then run on a separate node. E.g. on computerome:
qsub -W group_list=geogenetics -A geogenetics -l nodes=1:ppn=40,walltime=100:00:00,mem=100gb -N paint_panel0 -F 00 run_batch.sh

What this does: runs painting of one target individual per thread: 
bash paintsample1by1.sh UKBB0 1 4 ../cp_panel_scripts
[NB script paintsample1by1.sh must be edited if size of ref panel changes!]

This makes a new temporary directory for each individual, including phasefile with 1 target + all reference. It then runs will_04-paintvspanel.sh, which in turn runs will_paint_withinpanel-b.sh. This paints the target individuals and stores the local painting output in a memory efficient format. 

When these have all run, we cat the results together in the right order. Commands (where ukbb_samples is a file containing the 24000 IDs of the target individuals, one per line): 
python3
import pandas as pd
df=pd.read_csv("ordered_all_pop_ids_mapped", sep=" ", header=None)
df=df[0]
df1 = df.head(24000)
df1.to_csv("ukbb_samples", index=False, header=False, sep=" ")
exit()
chrlist=`seq 1 22`
for chr in $chrlist 
do while read p
do cat temp.$p/will_modernvsancient/painting/$chr.all_copyprobsperlocus.txt >> $chr.master_all_copyprobsperlocus.txt
done < ukbb_samples
done

while read p; do
awk "NR==3" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.chunkcounts.out >> ordered_all_pop_ids_mapped.allchr.chunkcounts.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.chunklengths.out >> ordered_all_pop_ids_mapped.allchr.chunklengths.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.mutationprobs.out >> ordered_all_pop_ids_mapped.allchr.mutationprobs.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out >> ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out >> ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out
done < ukbb_samples
mkdir perchrom_results
for chr in $chrlist
do while read p; do
awk "NR==3" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.chr$chr.chunkcounts.out >> perchrom_results/ordered_all_pop_ids_mapped.chr$chr.chunkcounts.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.chr$chr.chunklengths.out >> perchrom_results/ordered_all_pop_ids_mapped.chr$chr.chunklengths.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.chr$chr.mutationprobs.out >> perchrom_results/ordered_all_pop_ids_mapped.chr$chr.mutationprobs.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.chr$chr.regionchunkcounts.out >> perchrom_results/ordered_all_pop_ids_mapped.chr$chr.regionchunkcounts.out
awk "NR==2" temp.$p/will_modernvsancient/painting/ordered_all_pop_ids_mapped.chr$chr.regionsquaredchunkcounts.out >> perchrom_results/ordered_all_pop_ids_mapped.chr$chr.regionsquaredchunkcounts.out
done < ukbb_samples
done

Then run clean_and_repaint.sh which checks the results files, removes bad lines, and writes commands to repaint remaining individuals. E.g. on computerome:
qsub -W group_list=geogenetics -A geogenetics -l nodes=1:ppn=1,walltime=100:00:00,mem=100gb -N clean_and_repaint clean_and_repaint.sh

Once repainted (i.e. you’ve run the new commands in batch_files), run clean_and_repaint.sh again. Repeat until there are no more commands to run (this may happen on the first time, or it might take a few goes). 

General notes:
•	Output of local painting is reversed order compared to recombfiles (i.e. SNP position is decreasing). 
•	Dependencies: if you have fineSTRUCTURE installed and it is running without problem then that should be the main hurdle overcome. Check that you are using fs 4.0.1. The scripts use python3 with various mainstream packages (e.g. numpy, pandas) which shouldn’t cause issues. 
