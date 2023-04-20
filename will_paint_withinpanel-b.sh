#!/bin/bash
## paint_withinpanel.sh
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Copyright Daniel Lawson, 2020
## All rights reserved.

if [ "$#" -ne "3" ] ; then
    echo "Usage: paint_withinpanel.sh <full/repaint> <dir> <refdata.cp>"
    echo "<full/repaint>: whether we paint and repaint with a prior, or do a single painting"
    echo "<dir>: The root location of the output, which will be placed in <dir>/painting/<name>"
    echo "<refdata.cp>: A .cp file explaining the reference data: must contain the popidfile, donoridfile, recombfiles and phasefiles, popNeinf, popmuinf"
    exit 0
fi
set -e

myecho(){
    if [ $verbose ] ; then
	echo $@
    fi
}

verbose=TRUE
rmcp=TRUE
mode="$1"
dir="$2/painting"
refcp="$3"
cmdlist="$dir/commandlist1.txt"

#############
echo "Using directory $dir"
mkdir -p $dir
echo "Using cp file $refcp"
refpanelids=`grep ^popidfile $refcp | grep -o '^[^#]*' | cut -f2 -d:`
echo "Using reference ids file $refpanelids"
tmp=`echo $refpanelids | sed 's/.ids$//'`
refpanelname=`basename $tmp`
echo "Using reference panel name $refpanelname"
refpaneldonor=`grep ^donoridfile $refcp | grep -o '^[^#]*' | cut -f2 -d:`
myecho "Using reference donor file $refpaneldonor"
matchfile=`grep ^matchfile $refcp | grep -o '^[^#]*' | cut -f2 -d:`
if [ "$matchfile" != "" ] ;then
    myecho "Using matching file $matchfile"
else
    myecho "No matching file specified"
fi

grep " R$" $refpaneldonor | cut -f1 -d' ' > $dir/recipientpopulations.txt
cat $refpanelids | grep " 1$" | grep -Fwf $dir/recipientpopulations.txt > $dir/recipientindividuals.txt
numrecips=`cat $dir/recipientindividuals.txt | wc -l`
myecho "Counted $numrecips recipient individuals"

#############
cat $dir/recipientindividuals.txt | awk '{printf(" '$dir'/painting/%s/'$refpanelname'/%s.chunkcounts.out\t%s\n",$1,$1,$1);}' > $dir/allfiles.txt
touch $dir/testfileresults.txt
bash testcpoutput.sh $dir/allfiles.txt > $dir/testfileresults.txt

#############
if [ `cat $dir/testfileresults.txt | wc -l` -gt 0 ] ; then
    cat $dir/testfileresults.txt | awk '{print $2}' > $dir/idstoprocess.txt
    rm -f $cmdlist
    while read name; do
	if [ "$matchfile" != "" ] ;then
	    refname=`grep -w "^$name" $matchfile | cut -f2 -d' '`
	    if [ "$refname" == "" ] ; then
		refname="$name"
	    fi
	else
	    refname="$name"
	fi
	if [ "$refname" != "$name" ]; then
	    myecho "Using name $refname for name $name"
	fi
        echo "Test case $mode"
        echo bash will_paintsample_withinpanel-b.sh $mode $name $dir $refcp $refname >> $cmdlist
    done < $dir/idstoprocess.txt
    numruns=`cat $cmdlist | wc -l`
    echo "$numruns Runs required. Run with:"
    echo qsub_run.sh -f $cmdlist -n 8 -m 256
    exit 0
else
    echo "No problems detected, Running combine step"
fi
#############
nchromosomes=`grep ^phasefiles $refcp | grep -o '^[^#]*' | cut -f2 -d: | tr ',' '\n' | wc -l`
echo "Combine A: combine all whole genome results to $dir/$refpanelname.allchr"
fs combine -o $dir/$refpanelname.allchr `cat $dir/recipientindividuals.txt | awk '{printf(" '$dir'/painting/%s/'$refpanelname'/%s",$1,$1);}'` &> $dir/$refpanelname.allchr.log
if [ `tail -n 1 $dir/$refpanelname.allchr.log | grep "successfully" | wc -l` == "0" ]; then
    echo "ERROR: chromocombine error. See $dir/$refpanelname.allchr.log"
    exit 1
fi
## When this completes:
if [ "$mode" == "full" ] ;then
    ending="rerun"
else
    ending="$mode"
fi
for chr in `seq 1 $nchromosomes`; do
    echo "Combine B $chr of $nchromosomes: combine chromosome $chr results to $dir/$refpanelname.chr$chr"
    fs combine -o $dir/$refpanelname.chr$chr `cat $dir/recipientindividuals.txt | awk -v chr=$chr '{printf(" '$dir'/painting/%s/'$refpanelname'/cp/'$ending'.chr%s",$1,chr);}'` &> $dir/$refpanelname.chr$chr.log
    if [ `tail -n 1 $dir/$refpanelname.chr$chr.log | grep "successfully" | wc -l` == "0" ]; then
	echo "ERROR: chromocombine error. See $dir/$refpanelname.chr$chr.log"
	exit 1
    fi
done
echo "Success! After checking the chromocombine output, rm -r $dir/painting ."
