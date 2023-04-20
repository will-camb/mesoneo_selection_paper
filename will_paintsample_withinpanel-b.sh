#!/bin/bash
## paint_withinpanel.sh
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Copyright Daniel Lawson, 2020
## All rights reserved.

set -e
echo "Found $# arguments"
if [ "$#" -lt "4" ] ; then
    echo "Usage: will_paintsample_withinpanel-b.sh <full/repaint> <name> <dir> <refdata.cp> <optional: removename>"
    echo "<full/repaint>: whether we paint and repaint with a prior, or do a single painting"
    echo "<name>: The sample name to be processed"
    echo "<dir>: The root location of the output, which will be placed in <dir>/painting/<name>"
    echo "<refdata.cp>: A .cp file explaining the reference data: must contain the popidfile, donoridfile, recombfiles and phasefiles, popNeinf, popmuinf"
    echo "<removename>: An optional name of an individual in the panel that you want to treat as \"self\", and remove from the painting."
    exit 0
fi

myecho(){
    if [ $verbose ] ; then 
	echo $@
    fi
}

verbose=TRUE
rmcp=TRUE
mode="$1"
name="$2"
dir="$3/painting/$name"
refdatacp="$4"
if [ $# -gt 4 ]; then
    removename="$5"
else
    removename="$4"
fi

echo "Processing individual $name (Removing $removename)..."

echo "Reading $refdatacp"
refloc=`dirname $refdatacp`
myecho "Using reference location $refloc"
refpanelne=`grep ^popNeinf $refdatacp | grep -o '^[^#]*' | cut -f2 -d:`
refpanelmu=`grep ^popmuinf $refdatacp | grep -o '^[^#]*' | cut -f2 -d:`
myecho "Using reference mu $refpanelmu Ne $refpanelne"
refpanelids=`grep ^popidfile $refdatacp | grep -o '^[^#]*' | cut -f2 -d:`
myecho "Using reference ids file $refpanelids"
refpaneldonor=`grep ^donoridfile $refdatacp | grep -o '^[^#]*' | cut -f2 -d:`
myecho "Using reference donor file $refpaneldonor"
tmp=`echo $refpanelids | sed 's/.ids$//'`
refpanelname=`basename $tmp`
myecho "Using reference panel name $refpanelname"



if [ "$mode" == "full" ] ; then
    echo "Painting using uniform prior first, to obtain overall painting"
elif [ "$mode" == "repaint" ] ; then
    echo "Repainting using previous run as a prior"
elif [ "$mode" == "single" ] ; then
    echo "Painting using uniform prior only"
else
    echo "ERROR: Unknown mode $mode"
    exit 1
fi

mkdir -p "$dir/$refpanelname/cp/"

## Extract data
echo "$name $name 1" > "$dir/$refpanelname/cp/sample.ids"

nchromosomes=`grep ^phasefiles $refdatacp | grep -o '^[^#]*' | cut -f2 -d: | tr ',' '\n' | wc -l`
echo "Counted $nchromosomes chromosomes" 
chrlist=`seq 1 $nchromosomes`
if [ "$mode" == "full" ] || [ "$mode" == "single" ] ; then
    output="initial"
    for chr in $chrlist ; do
	echo "Processing Chr $chr"
	#Rscript makeIds.R -n "$name" "$refpanelids" "$dir/$refpanelname/initial.chr$chr.ids" "$removename"
	Rscript makeIds.R -k $refpanelids $dir/$refpanelname/initial.chr$chr.ids $removename
	echo "Done makeIds.R"
	Rscript orderRecipientsByDonorFile.R $dir/$refpanelname/initial.chr$chr.ids $refpaneldonor $dir/$refpanelname/donororder.initial.chr$chr.ids
	echo "Done orderRecipientsByDonorFile.R"
	echo "Using name $name"
	echo getindividualnumverfromname.sh $name $dir/$refpanelname/donororder.initial.chr$chr.ids 1 $refpaneldonor
	indnum=`bash getindividualnumverfromname.sh $name $dir/$refpanelname/donororder.initial.chr$chr.ids 1 $refpaneldonor`

	if [ "$indnum" == "" ] ;then
	    echo "ERROR!"
	    exit 1
	fi
	echo "Found individual named $name at number $indnum"
	refpanelphase=`grep ^phasefiles $refdatacp | grep -o '^[^#]*' | cut -f2 -d: | cut -f$chr -d','`
	refpanelrecomb=`grep ^recombfiles $refdatacp | grep -o '^[^#]*' | cut -f2 -d: | cut -f$chr -d','`
## Extract data
	myecho "ReferencePhase: $refpanelphase ReferenceRec: $refpanelrecomb IndNum: $indnum"
	cmd="fs cp -n $refpanelne -M $refpanelmu -g $refpanelphase -r $refpanelrecomb -f $refpaneldonor $indnum $indnum -t $dir/$refpanelname/initial.chr$chr.ids -o $dir/$refpanelname/cp/$output.chr$chr -s 0"
	myecho $cmd
	$cmd
	echo "Done extract data phase, now making the prior"
    done
## Make the prior, i.e. the copying averages
    if [ $mode == "full" ] ; then
	Rscript makePrior.R -o $dir/$refpanelname/cp/$output $dir/$refpanelname/cp/$output.chr*.chunklengths.out
	echo "Done makePrior.R"
    fi
fi

if [ "$mode" == "full" ] || [ "$mode" == "repaint" ] ;then
    output="rerun"
    for chr in $chrlist ; do
	outdat="$dir/$refpanelname/complete.chr$chr.phase"
	echo "Processing Chr $chr"
	Rscript makeIds.R -k $refpanelids $dir/$refpanelname/rerun.chr$chr.ids $removename
	Rscript orderRecipientsByDonorFile.R $dir/$refpanelname/rerun.chr$chr.ids $refpaneldonor $dir/$refpanelname/donororder.rerun.chr$chr.ids
	indnum=`bash getindividualnumverfromname.sh $name $dir/$refpanelname/donororder.rerun.chr$chr.ids 1 $refpaneldonor`
	if [ "$indnum" == "" ] ;then
	    echo "ERROR!"
	    exit 1
	fi
	echo "Found individual named $name at number $indnum"
	refpanelphase=`grep ^phasefiles $refdatacp | grep -o '^[^#]*' | cut -f2 -d: | cut -f$chr -d','`
	refpanelrecomb=`grep ^recombfiles $refdatacp | grep -o '^[^#]*' | cut -f2 -d: | cut -f$chr -d','`
## Extract data
	echo "Making initial all_copyprobsperlocus.out at $3 to add each copyprobsperlocus.out to as they're created"
  python make_all_copyprobsperlocus.out.py -chr $chr -phasefile $refpanelphase -o $3
	myecho "ReferencePhase: $refpanelphase ReferenceRec: $refpanelrecomb IndNum: $indnum"
	cmd="fs cp -b -n $refpanelne -M $refpanelmu -g $refpanelphase -r $refpanelrecomb -f $refpaneldonor $indnum $indnum -t $dir/$refpanelname/rerun.chr$chr.ids -o $dir/$refpanelname/cp/rerun.chr$chr -s 0"
	echo "Reprocessing Chr $chr"
	myecho $cmd
	$cmd
	gunzip $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out.gz
  python modify_copyprobsperlocus.out.py -copyprobsperlocus_location $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out
  cat $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out | sed 's/ //g' > $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out_modified
	cat $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out_modified >> $3/$chr.all_copyprobsperlocus.txt
	rm -r $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out $dir/$refpanelname/cp/rerun.chr$chr.copyprobsperlocus.out_modified
    done
fi

if [ "$nchromosomes" == "1" ] ;then
    fs combine -o $dir/$refpanelname/$name $dir/$refpanelname/cp/$output.chr1
else
    fcmd=`eval echo $dir/$refpanelname/cp/$output.chr{1..$nchromosomes}.chunkcounts.out`
    files=`ls $fcmd`
    fs combine -o $dir/$refpanelname/$name $files
fi

## Tidy up the temporary files
#rm -r $dir/$refpanelname/cp/
