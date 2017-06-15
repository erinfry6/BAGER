## read in tissue argument

tissue<-commandArgs()

## set paths to directories, be sure to modify your home directory and the Anccestral Reconstruction directory you are analyzing
path="/Users/lynchlab/Desktop/ErinFry/workflowr/AGER/"
pathData=paste(path,"data/forBAGER/",sep="")
pathResults=paste(path,"data/BAGERresults/",tissue,"/",sep="")
pathSSS=paste(pathResults,"SSS/",sep="")
pathKappa=paste(pathSSS,"kappa/",sep="")
pathDelta=paste(pathSSS,"delta/",sep="")
pathNone=paste(pathSSS,"none/",sep="")
pathKD=paste(pathSSS,"kd/",sep="")

######################################################################
#### define the genes to find the best model for
setwd(pathKappa)
options(stringsAsFactors = FALSE)
ldf <- list() # creates a list
listcsv<-paste("gene",as.character((order(dir(pattern = "*.txt")))),".txt", sep="")

## set the number of genes to evaluate
numgenes<-length(listcsv)

#collect the log marginal likelihood of the stepping stone sampler for each model for each gene using the following function: 

findLML<-function(numgenes){
  options(stringsAsFactors = FALSE)

  #find the log marginal likelihood for each gene
  logmarginallikelihood<-vector(length=length(ldf))
  for (i in 1:numgenes){
    logmarginallikelihood[i]<-read.csv(listcsv[i], sep='\t')[nrow(read.csv(listcsv[i], sep='\t')),]
  }
  
  return(logmarginallikelihood)
  
}


## find LML under each model
#### KAPPA
setwd(pathKappa)
KappaLML<-findLML(numgenes)


#### DELTA
setwd(pathDelta)
DeltaLML<-findLML(numgenes)

#### None
setwd(pathNone)
NoneLML<-findLML(numgenes)

#### KAPPA and DELTA
setwd(pathKD)
KDLML<-findLML(numgenes)


######################################################################


## Identify the best model for each gene
## choice of 1=Kappa, 2= KD, 3=Delta, 4=none
choice=vector()
for (i in 1:numgenes){
choice[i]<-which.max(c(KappaLML[i],KDLML[i],DeltaLML[i], NoneLML[i]))
}


## import the gene information to we can include them in our analysis
setwd(pathResults)
write.table(t(choice),"modelchoice.txt",col.names = FALSE, row.names = FALSE)

######################################################################


