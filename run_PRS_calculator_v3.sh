#!/bin/bash
#cd $PBS_O_WORKDIR
#module load anaconda3/4.4.0
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "5" ] ; then
    echo "Usage: run_PRS_calculator.sh <phenotype_file> <copyprobs_directory> <phasefile_directory> <idfile>"
    echo "<phenotype_file>: file containing filenames of phenotype GWAS results, 1 per row eg E4_OBESITY.gwas.imputed_v3.both_sexes.tsv.bgz"
    echo "<copyprobs_directory>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1')"
    echo "<bootstrap True/False>: Whether to bootstrap or not. Can take value either True or False"
    exit 0
fi

phenotype_file="$1"
copyprobs_directory="$2"
phasefile_directory="$3"
idfile="$4"
bootstrap="$5"
chrlist=$(seq 1 22)
ancestries="CHG EHG Farmer African EastAsian WHG Yamnaya"

# For running on a subset of painted individuals in UKB
for chr in $chrlist; do
  for anc in $ancestries; do
    echo "python3 PRS_calculator_v3.py -copyprobs_file $copyprobs_directory/temp.$anc.$chr.master_all_copyprobsperlocus.txt.gz -phasefile $phasefile_directory/transformed.$chr.merged.phase.gz -idfile $idfile -phenotypes $phenotype_file -chr $chr -anc $anc -bootstrap $bootstrap" >> PRS_calculator_commands_temp
  done
done

# For running on all painted individuals in UKB
#for chr in $chrlist; do
#  for anc in $ancestries; do
#    echo "python3 PRS_calculator_v3.py -copyprobs_file $copyprobs_directory/$anc.$chr.master_all_copyprobsperlocus.txt.gz -phasefile $phasefile_directory/$chr.merged.phase.gz -idfile $idfile -phenotypes $phenotype_file -chr $chr -anc $anc -bootstrap $bootstrap" >> PRS_calculator_commands_temp
#  done
#done


echo "Now running commands in PRS_calculator_commands in parallel"
shuf PRS_calculator_commands_temp > PRS_calculator_commands_"$phenotype_file"
rm PRS_calculator_commands_temp
cat PRS_calculator_commands_"$phenotype_file" | parallel  -j 50 # Change as appropriate

echo "*** All done, results are in PRS_calculations_v3 ***"

rm PRS_calculator_commands_"$phenotype_file"