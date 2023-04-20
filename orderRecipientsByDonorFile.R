#!/opt/anaconda3/bin/Rscript
## orderRecipientsByDonorFile.R
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Copyright Daniel Lawson, 2020
## Released under GPLv3
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   http://www.gnu.org/licenses

options(warn=-1)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<3){
    cat("Usage: Rscript orderRecipientsByDonorFile.R <input id file> <input donor file> <output id file>
Generate an ID file ordered by the population recipient labels in the order given in the input donor file.
")
    stop("Must provide 3 arguments")
}

tempargs=args

idf=tempargs[1]
df=tempargs[2]
idfout=tempargs[3]

idtab<-read.table(idf,as.is=TRUE,na.strings="N/A")
donortab<-read.table(df,as.is=TRUE,na.strings="N/A")
donors=donortab[donortab[,2]=="R", 1]
names(donors)=donors

idtab=idtab[idtab[,3]==1,,drop=FALSE]

idtablist=lapply(donors,function(x){
  idtab[idtab[,2]==x,,drop=FALSE]
})

idtabnew=do.call(rbind,idtablist)

write.table(idtabnew,file=idfout,quote=FALSE,
            row.names=FALSE,col.names=FALSE)
