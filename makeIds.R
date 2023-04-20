#!/opt/anaconda3/bin/Rscript
options(warn=-1)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<2){
    cat("Usage: Rscript makeIds.R <options> <input id file> <output id file> <optional: indname>
Generate a file with 1 sample removed from each population.
options:
-s : seed
-k : keep self in the population (for donor-mode use in finestructure)
<indname>: make sure we remove this individual from the population
")
    stop("Must provide 2 arguments")
}

tempargs=args
verbose=FALSE
seed=0
rmself=TRUE

i=1
while(i<length(tempargs)){
  if(tempargs[i] == "-s"){
    seed=as.numeric(tempargs[i+1] )
    tempargs=tempargs[- (i+(0:1)) ]
    }else if(tempargs[i]=="-k"){
      tempargs=tempargs[- i ]
      rmself=FALSE
    }else{
      i=i+1
    }
}
  
if(seed>0) set.seed(seed)

idf=tempargs[1]
idfout=tempargs[2]
id=""
if(length(tempargs)>2) id=tempargs[3]

idtab<-read.table(idf,as.is=TRUE,na.strings="N/A")

tind=which(idtab[,1]==id)
tpop=""
if( length(tind)>0 && !is.na(idtab[tind,2]) ) tpop=idtab[tind,2]

allpops=unique(idtab[idtab[,3]==1,2])

for (pop in allpops){
#  print(paste("pop=",pop,"tpop=",tpop))
  if(tpop==pop){
    if(rmself){
      myrem=tind
    }else{
      myrem=numeric()
    }
  }else{
    tw=which(idtab[,2]==pop & idtab[,3]==1)
    if(length(tw)==1) stop("Can't have only 1 individual in a population")
    myrem=sample(tw,1)
  }
    idtab[myrem,3]=0
}

write.table(idtab,file=idfout,quote=FALSE,
            row.names=FALSE,col.names=FALSE)

