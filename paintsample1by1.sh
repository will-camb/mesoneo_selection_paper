#!/bin/bash

echo "Found $# arguments"
if [ "$#" -lt "4" ] ; then
    echo "Usage: paintsample1by1.sh <id> <n> <phaselinenumber> <cp_panel_scripts>"
    echo "<id>: the individual which is being painted against the reference panel"
    echo "<n>: the number of the individual in the idfile"
    echo "<phaselinenumber>: the line of the phasefile containing the individual's first haplotype"
    echo "<cp_panel_scripts>: the location of the panel scripts which will be copied to each temp directory"
    exit 0
fi
set -e

name="$1"
number="$2"
phaselinenumber="$3"
cp_panel_scripts="$4"
dir="temp.$name"
chrlist=`seq 1 22`
nhaps=638
#module load tools finestructure/4.1.1
#module load parallel/20200922
#module load anaconda3/4.4.0
#module load intel/perflibs
#module load gcc
#module load R/3.6.1

for chr in $chrlist; do
  if [ ! -f $chr.master_all_copyprobsperlocus.txt ]; then
    echo "Making $chr.master_all_copyprobsperlocus.txt"
    touch $chr.master_all_copyprobsperlocus.txt
  fi
  done

if [ ! -f ordered_all_pop_ids_mapped.allchr.chunkcounts.out ]; then
    echo "Making ordered_all_pop_ids_mapped.allchr.chunkcounts.out, ordered_all_pop_ids_mapped.allchr.chunklengths.out, ordered_all_pop_ids_mapped.allchr.mutationprobs.out, ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out, ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out"
    touch ordered_all_pop_ids_mapped.allchr.chunkcounts.out
  fi
if [ ! -f ordered_all_pop_ids_mapped.allchr.chunklengths.out ]; then
    touch ordered_all_pop_ids_mapped.allchr.chunklengths.out
  fi
if [ ! -f ordered_all_pop_ids_mapped.allchr.mutationprobs.out ]; then
    touch ordered_all_pop_ids_mapped.allchr.mutationprobs.out
  fi
if [ ! -f ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out ]; then
    touch ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out
  fi
if [ ! -f ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out ]; then
    touch ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out
  fi

if [ -d $dir ]; then
    echo "Aborting this individual as $dir already exists"
    exit 1
  fi

echo "Using temporary directory $dir"
mkdir -p "$dir"
mkdir -p "$dir/phasefiles"
cmd="cp $cp_panel_scripts/* $dir"
echo "Copying scripts from $cp_panel_scripts"
$cmd
echo "Done copying scripts, now making new idfile at $dir"
{ head -n $number ordered_all_pop_ids_mapped | tail -n 1 && tail -n 318 ordered_all_pop_ids_mapped; } > $dir/ordered_all_pop_ids_mapped

phaselinenumber2=$(($phaselinenumber + 1))
for chr in $chrlist ; do
  touch $dir/phasefiles/$chr.merged.phase
  head -n 3 phasefiles/$chr.merged.phase > $dir/phasefiles/$chr.merged.phase
  awk "NR>=$phaselinenumber && NR<=$phaselinenumber2" phasefiles/$chr.merged.phase >> $dir/phasefiles/$chr.merged.phase
  tail -n 636 phasefiles/$chr.merged.phase >> $dir/phasefiles/$chr.merged.phase
#  sed -i '' "1s/.*/$nhaps/" $dir/phasefiles/$chr.merged.phase
  sed -i "1s/.*/$nhaps/" $dir/phasefiles/$chr.merged.phase
  echo "Copied lines to $dir/phasefiles/$chr.merged.phase"
  done

echo "Moving to directory $dir to run remaining commands"
cd $dir
cmd="bash will_04-paintvspanel.sh"
echo $cmd
$cmd
cmd="bash will_modernvsancient/painting/commandlist1.txt"
echo $cmd
$cmd
bash will_04-paintvspanel.sh
cd ../

echo "Now copying all_copyprobsperlocus.txt and chunkcounts/chunklengths etc to master files"
for chr in $chrlist ; do
  cat $dir/will_modernvsancient/painting/$chr.all_copyprobsperlocus.txt >> $chr.master_all_copyprobsperlocus.txt
  done

awk "NR==3" $dir/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.chunkcounts.out >> ordered_all_pop_ids_mapped.allchr.chunkcounts.out
awk "NR==2" $dir/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.chunklengths.out >> ordered_all_pop_ids_mapped.allchr.chunklengths.out
awk "NR==2" $dir/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.mutationprobs.out >> ordered_all_pop_ids_mapped.allchr.mutationprobs.out
awk "NR==2" $dir/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out >> ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out
awk "NR==2" $dir/will_modernvsancient/painting/ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out >> ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out

echo "Deleting temp $dir"
rm -r $dir