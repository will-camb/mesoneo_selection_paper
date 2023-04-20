#!/bin/bash
## getindividualnumverfromname.sh
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Copyright Daniel Lawson, 2020
## Released under GPLv3
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   http://www.gnu.org/licenses

help (){
    echo "Usage: getindividualnumverfromname.sh <name> <refidfile> <mode> <optional:donorfile>"
    echo "         <name> Give the name for the desired individual, i.e. the first field in the ID file, from which an individual number is extracted and used as Mode 1."
    echo "         <refidfile> Give the filename of the ID file for which we seek the individual number"
    echo "         <mode> should be either 0 or 1. If mode==0, we count ALL LINES in the file. If mode==1, we only count those lines that are included,i.e. ther 3rd field is 1."
    echo "         <donorfile> donor file listing RECIPIENT populations that to be counted. (Donors may be present but are ignored.)"
    echo "Outputs the number for the individual. Mode=1 should be used for use with a donorfile."
}


if [ "$#" -lt "3" ] ; then
    help
    exit 0
fi
indname="$1"
refidfile="$2"
mode="$3"
donorfile=""
if [ "$#" -eq "4" ] ; then
    donorfile="$4"
fi

if [ "$mode" == "0" ]; then
    if [ "$donorfile" == "" ]; then
	i=`grep -n "^$indname " $refidfile | cut -d: -f1`
    else
	(>&2 echo "ERROR: No usage case for mode=0 and donorfile specified! Not implemented")
	exit 1;
    fi
else
    if [ "$donorfile" == "" ]; then
	i=`cat $refidfile | grep " 1$" | grep -n "^$indname " | cut -d: -f1`
    else
	grep " R$" $donorfile | cut -f1 -d' ' > $indname.donorfile.tmp
	i=`cat $refidfile | grep " 1$" | grep -Fwf $indname.donorfile.tmp | grep -n "^$indname " | cut -d: -f1`
	rm $indname.donorfile.tmp
	if [ "$i" == "" ] ; then
	    (>&2 echo "ERROR: Individual $indname not found in ID file $refidfile as a recipient from file $donorfile.")
	    exit 1;
	fi
    fi
fi

if [ "$i" == "" ] ; then
    (>&2 echo "ERROR: Individual $indname not found in ID file $refidfile." ) 
    exit 1;
fi
echo $i