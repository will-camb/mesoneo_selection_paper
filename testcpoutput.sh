#!/bin/sh
## testcpoutput.sh
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Copyright Daniel Lawson, 2020
## Released under GPLv3
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   http://www.gnu.org/licenses

if [ $# -eq 0 ]; then
    echo "Usage: testcpoutput.sh <filename>"
    echo "Checks every line of the file for whether it is plausible CP output. If not, it outputs the file to stdout."
    echo "file can contain mutliple tab separated values; only the first is used. The whole line is written to file."
    exit 1
fi
infile=$1
while read p; do
    pp=`echo $p | awk '{print $1}'`
    if [ -e "$pp" ]; then
	nlines=`cat $pp | wc -l`
	if [ "$nlines" -lt 2 ]; then
	    echo $p
	fi
    else
	echo $p
    fi
done < $infile
