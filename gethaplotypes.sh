#!/bin/bash
help (){
    echo "Usage: gethaplotypes.sh <options>"
    echo "Has two modes: 1. use -i and specify the individual number. 2. use -r and -n and specify an external id file and the name."
    echo "Options:"
    echo "         -O <output.phase>: REQUIRED. Specify the output file location."
    echo "         -f <input.phase> : REQUIRED. Set the input filename."
    echo "         -i <i> : MODE 1 REQUIRED. extract the i'th individual."
    echo "         -r <reference id file> MODE 2 REQUIRED. Specify an ID file for the phase file."
    echo "         -n <name> MODE 2 REQUIRED. Give the name for the desired individual, i.e. the first field in the ID file, from which an individual number is extracted and used as Mode 1."
    echo "         -p <ploidy> (default: 2) Set the number of haplotypes per individual (remember we provide individual number as input)"
    echo "         -o : output haplotypes only, ie no header line or SNPS line"
    echo "         -h : show this help"
    echo "Outputs the contents of a valid phase file including header line, from <input.phase>" 
}

ploidy=2
output="full"
i=1
infile=""
outfile=""
refidfile=""
indname=""
while [ $# -gt 0 ]; do
    case $1 in
	-O) outfile=$2; shift;shift;;
	-i) i=$2; shift;shift;;
	-f) infile=$2; shift;shift;;
	-r) refidfile=$2; shift;shift;;
	-n) indname=$2; shift;shift;;
	-p) ploidy="$2"; shift;shift;;
	-o)  output="haplotypes"; shift;;
	-h) help;
	    exit 0;;
	*) help; exit 1;;
    esac
done

## Check if we've been given an input file
if [ "$infile" == "" ]; then
    help; exit 1
fi

## Check if we've been given an output file
if [ "$outfile" == "" ]; then
    help; exit 1
fi

## Check if we've been given an ID file
if [ "$refidfile" != "" ] ; then
    if [ "$indname" != "" ] ; then
	## Get the individual number from the ID file
	i=`getindividualnumverfromname.sh $indname $refidfile 0`
	if [ "$i" == "" ] ; then
	    echo "ERROR: Individual $indname not found in ID file $refidfile."
	    exit 1;
	fi
    else
	help; exit 1;
    fi
fi
##echo "ploidy=$ploidy,output=$full;infile=$infile;"
if [ "$i" == "" ] ; then
    echo "ERROR: No individual number specified or read from ID file. Use either -i or -r and -n as specified in the help."
    exit 1;
fi

haps=`echo "3+$i*$ploidy-$ploidy+1"| bc`
hapsmax=`echo "$haps + $ploidy - 1" | bc`
hapsl=`seq $haps $hapsmax`
hapsl=`echo $hapsl | paste -s | tr ' ' ','`
#hapsl=`echo $hapsl | tr '\n ' ','`

rm -f $outfile
touch $outfile
if [ "$output"=="full" ] ; then
    echo 2 >> $outfile
    sed -n 2,3p $infile >> $outfile
fi
sed -n ${hapsl}p $infile >> $outfile
