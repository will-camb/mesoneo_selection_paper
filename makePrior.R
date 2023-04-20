#!/cm/shared/languages/R-3.0.2/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<2){
      cat("Usage: Rscript makePrior.R <options> -o <output root> <list of single roots>
options:
")
          stop("Must provide 2 arguments")
    }

tempargs=args
outroot=""
if(tempargs[1] == "-o"){
  outroot=tempargs[2]
  tempargs=tempargs[- (1:2) ]
}
print(paste("Using output ",outroot))
if(outroot=="") stop("Invalid output root provided")
dirout=dirname(outroot)

allsingleroots=tempargs
allroots=paste(allsingleroots,collapse=" ")
system(paste0("mkdir -p ",dirout))

system(paste0("fs combine -o ",outroot," ",allroots))

tdata=read.table(paste0(outroot,".chunklengths.out"),header=T,row.names=1)
odata<-data.frame(names=c(colnames(tdata)),dr=c(rep("D",dim(tdata)[2])),prior=c(as.numeric(tdata/sum(tdata))))
write.table(odata,file=paste0(dirout,"/prior.donor"),quote=FALSE,col.names=FALSE,row.names=FALSE)
