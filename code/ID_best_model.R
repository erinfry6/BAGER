tissue="br"

## set paths to directories, be sure to modify your home directory and the Anccestral Reconstruction directory you are analyzing
path="/Users/lynchlab/Desktop/ErinFry/workflowr/AGER/"
pathData=paste(path,"data/forBAGER/",sep="")
pathResults=paste(path,"data/BAGERresults/",tissue,"/",sep="")
pathSSS=paste(pathResults,"Stones/",sep="")

## set the models tested
models<-c("VarRates", "Lambda")

######################################################################

## read in the list of genes to see which model is best
options(stringsAsFactors = FALSE)
setwd(paste(pathSSS, models[length(models)],sep=""))
listcsv<-paste("gene",as.character((order(dir(pattern = "*.txt")))),".txt", sep="")

## create matrix to hold log marginal information for each gene and model
LML<-matrix(nrow = (length(models)+2), ncol = length(listcsv))
row.names(LML)<-c("gene.number",models,"model.choice")
LML[1,]<-listcsv[1:ncol(LML)]

# for each model and gene, pull out the log marginal likelihood
for (m in 1:length(models)){
  setwd(paste(pathSSS, models[m],sep=""))
  LML[1+m,]<-sapply(listcsv, function(x){read.csv(x, sep='\t')[nrow(read.csv(x, sep='\t')),]})
}

## for each gene, find the best model
for (i in 1:ncol(LML)){
  LML[nrow(LML),i]<-models[which.max(LML[2:(nrow(LML)-1),i])]
}

## save this information in tissue results directory
write.table(LML,paste(pathResults,"modelchoice.txt", sep=""),col.names = FALSE, sep='\t')

######################################################################



