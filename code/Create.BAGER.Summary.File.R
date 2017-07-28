
## Erin Fry
## July 28 2017

## this script identifies the best model for each gene's evolution

## MODIFICATION INSTRUCTIONS:
## all required or recommended changes are in the first section of the script
## change the variable tissue to the tissue code you are evaluating
## change the variable path to home directory for the project
## if testing other models, modify model vector to include all models
## you must modify the nodes of interest 'if then' section. Change the names of the nodes to the clade name 
## and the "Node ..." to the column names in the files created by BayesTraits for the .tree for that tissue

######################################################################

## SET PATHS, ARGUMENTS AND LIBRARIES

tissue<-"br"

## set paths to directories, be sure to modify your home directory and the Anccestral Reconstruction directory you are analyzing
path="/Users/lynchlab/Desktop/ErinFry/workflowr/AGER/"
pathData=paste(path,"data/forBAGER/",sep="")
pathResults=paste(path,"data/BAGERresults/",tissue,"/",sep="")
pathAncRecon=paste(pathResults,"AncStates/",sep="")

## set the models tested
models<-c("VarRates", "Lambda")

## set the nodes of interest
if (tissue=="br"){
  node.ancPrimates="Node-00000"
  node.ancHaplorhini="Node-00001"
  node.ancApes="Node-00002"
  node.ancHominini="Node-00003"
  node.ancHomo="Node-00004"
  node.ancPan="Node-00008"
  node.ancPtr="Node-00009"
  node.ancPpa="Node-00016"
  node.ancGorilla="Node-00020"
  node.ancOrang="Node-00023"
  node.ancMacaque="Node-00026"
  
  
} else if (tissue=="cb"){
  node.ancPrimates="Node-00000"
  node.ancHaplorhini="Node-00001"
  node.ancApes="Node-00002"
  node.ancHominini="Node-00003"
  node.ancHomo="Node-00004"
  node.ancPan="Node-00008"
  node.ancPtr="Node-00009"
  node.ancPpa="Node-00016"
  node.ancGorilla="Node-00020"
  node.ancOrang="Node-00023"
  node.ancMacaque="Node-00026"
  
} else if (tissue=="lv"){
  node.ancPrimates="Node-00000"
  node.ancHaplorhini="Node-00001"
  node.ancApes="Node-00002"
  node.ancHominini="Node-00003"
  node.ancHomo="Node-00004"
  node.ancPan="Node-00008"
  node.ancPtr="Node-00009"
  node.ancPpa="Node-00016"
  node.ancGorilla="Node-00020"
  node.ancOrang="Node-00023"
  node.ancMacaque="Node-00026"
  
} else if (tissue=="kd"){
  node.ancPrimates="Node-00000"
  node.ancHaplorhini="Node-00001"
  node.ancApes="Node-00002"
  node.ancHominini="Node-00003"
  node.ancHomo="Node-00004"
  node.ancPan="Node-00008"
  node.ancPtr="Node-00009"
  node.ancPpa="Node-00016"
  node.ancGorilla="Node-00020"
  node.ancOrang="Node-00023"
  node.ancMacaque="Node-00026"
  
} else if (tissue=="ts"){
  node.ancPrimates="Node-00000"
  node.ancHaplorhini="Node-00001"
  node.ancApes="Node-00002"
  node.ancHominini="Node-00003"
  node.ancHomo="Node-00004"
  node.ancPan="Node-00008"
  node.ancPtr="Node-00009"
  node.ancPpa="Node-00016"
  node.ancGorilla="Node-00020"
  node.ancOrang="Node-00023"
  node.ancMacaque="Node-00026"
  
} else if (tissue=="ht"){
  node.ancPrimates="Node-00000"
  node.ancHaplorhini="Node-00001"
  node.ancApes="Node-00002"
  node.ancHominini="Node-00003"
  node.ancHomo="Node-00004"
  node.ancPan="Node-00008"
  node.ancPtr="Node-00009"
  node.ancPpa="Node-00016"
  node.ancGorilla="Node-00020"
  node.ancOrang="Node-00023"
  node.ancMacaque="Node-00026"
  
} else {
  print("this tissue is not known") }

## load libraries
library(dplyr)
library(biomartr)

## set strings as factors to false so can read in properly
options(stringsAsFactors = FALSE)

######################################################################
## DEFINE FUNCTIONS

## define function for reading in these types of files which will read incorrectly using standard read.csv
read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

## Import the ancestral expression reconstructions, excluding the header created by BayesTraits run with '12 VarRates'
read.AncRecon=function(file, firstrow, lastrow,header=F, sep='\t'){
  temp<-t(read.tcsv(file, sep='\t', header=F)[,firstrow:lastrow])
  colnames(temp)<-temp[1,]
  temp<-temp[-1,]
  return(temp)
}

## Calculate the Bayesian Posterior Probability of Divergnce (BPPD)
## first, find the difference between the two reconstructions
## then, take the maximum proportion of the iterations that are greater than or less than 0
calc.BPPD=function(file, recon1, recon2){
  diff<-as.numeric(file[,which(colnames(file)==paste(recon1, " - 1",sep=""))])-as.numeric(file[,which(colnames(gene)==paste(recon2, " - 1",sep=""))])
  BPPD<-max(1-(length(which(diff>0))/(expectedrows-1)), (length(which(diff>0))/(expectedrows-1)))  ## to do
  return(BPPD)
}


## COLLECT STATS FUNCTION

collect.stats.for.one.gene<-function(gene, reconAnc, reconDesc){
  
  ## calculate the fold difference increase in the ancestral node compared to the descendent node
  foldSD<-sd(gene[,which(colnames(gene)==paste(reconAnc, " - 1",sep=""))])/sd(gene[,which(colnames(gene)==paste(reconDesc, " - 1",sep=""))])
  
  ## calculate bppd
  BPPD<-calc.BPPD(gene, reconAnc, reconDesc)
  
  ## find the median reconstructed values for both nodes
  MedianAnc<-median(as.numeric(gene[,which(colnames(gene)==paste(reconAnc, " - 1",sep=""))]))
  MedianDesc<-median(as.numeric(gene[,which(colnames(gene)==paste(reconDesc, " - 1",sep=""))]))
  
  ## return all stats for that gene and the pair of ancestral reconstructions
  return(c(foldSD,BPPD,MedianAnc,MedianDesc))
}


######################################################################
## IMPORT GENE DATA
## import the gene information to we can include them in our analysis
genenames<-read.csv(paste(pathData,tissue,"_genesincluded.txt",sep=""),header=T, sep='\t')

## import other infomration of those genes from biomart
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(genenames)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol","ensembl_gene_id","chromosome_name"),values=genes,mart= mart)

## because the first gene is duplicated, duplicate the first gene in G_list to stay consistant
G_list<-rbind(G_list[1,],G_list)

## load model choice information
modelchoice<-t(read.csv(paste(pathResults,"modelchoice.txt", sep=""),header=F, sep="\t"))
colnames(modelchoice)<-modelchoice[1,]
modelchoice<-as.data.frame(modelchoice[-1,])

######################################################################
## SAVE DETAILS ABOUT BAYESTRAITS RESULTS FILES

## define the genes to be summarized
listcsv<-modelchoice$gene.number

## sometimes there are multiple ensembl id's per gene name, and this will misalign gene labels with the reconstruction information
dup.ensembl<-which(rownames(genenames) %in% G_list$ensembl_gene_id ==FALSE)

if (length(dup.ensembl)>0 ){
  warning(paste(length(dup.ensembl), "genes are duplicated v ensembl"))
  ## eliminate those genes from the list (it should be very few genes)
  listgenes<-listcsv[-dup.ensembl]
} else {
  listgenes<-listcsv
}


## find which row the iteration information begins on for this tissue's tree
finding.information.about.file<-(read.tcsv(paste(pathAncRecon,models[length(models)],"/",listcsv[1],sep=""), sep='\t'))
it.begin<-which(colnames(finding.information.about.file)=="Itter")

## set the total number of rows you expect to have in each file so the code will warn you if reconstructions failed
expectedrows=ncol(finding.information.about.file)-it.begin


######################################################################
## PREPARE SUMMARY DATA FRAME
## create the summary dataframe, first with gene information
summary.df<-as.data.frame(cbind(listgenes, modelchoice$model.choice[1:length(listgenes)], G_list[1:length(listgenes),]))
colnames(summary.df)<-c("listcsv","modelchoice","hgnc_symbol", "ensembl_gene_id", "chromosome_name")

## define the statistics to be collected 
stats.to.collect<-c("foldSD","BPPD","reconAnc","reconDesc")

## and each lineage in the tree
lineages.to.test<-matrix(ncol=2, byrow = TRUE,data=c("ancPrimates", "ancMacaque",
                                                     "ancPrimates", "ancHaplorhini",
                                                     "ancHaplorhini", "ancOrang",
                                                     "ancHaplorhini","ancApes",
                                                     "ancApes","ancGorilla",
                                                     "ancApes", "ancHominini",
                                                     "ancHominini","ancPan",
                                                     "ancHominini", "ancHomo",
                                                     "ancPan", "ancPtr",
                                                     "ancPan", "ancPpa"))

## define column numbers for completing the dataframe with statistics
basecol<-ncol(summary.df)
statcol<-length(stats.to.collect)

## create the columns for each lineage's statistics
for (l in 1:nrow(lineages.to.test)){
  for (x in stats.to.collect){
    newcolname<-paste(lineages.to.test[l,1],lineages.to.test[l,2],x,sep="_")
    summary.df[[newcolname]]<-"na"
  }
}

######################################################################
## COMPLETE THE SUMMARY DATA FRAME BY COLLECTING STATISTICS FOR EACH GENE AND EACH LINEAGE
## for each gene
for (i in 1:length(listgenes)){
  
  ## load the model choice
  choice<-modelchoice$model.choice[i]
  
  ## import reconstruction information
  gene<-read.AncRecon(paste(pathAncRecon,choice,"/",listcsv[i],sep=""), firstrow = it.begin,lastrow=(it.begin+expectedrows), sep='\t') # read in the output of BayesTraits
  if (nrow(gene)!=expectedrows){
    warning(paste("The reconstruction of gene number", i, "failed. Check file"))
  }
  
  ## for each lineage
  for (l in 1:nrow(lineages.to.test)){
    ## find the BAGER statistics
    summary.df[i,(basecol+(l-1)*statcol+(1:statcol))]<-collect.stats.for.one.gene(gene, reconAnc = eval(as.symbol(paste("node.",lineages.to.test[l,1],sep=""))), eval(as.symbol(paste("node.",lineages.to.test[l,2],sep=""))))
    
  }


}

## view summary stats
head(summary.df)

######################################################################
## SAVE BAGER SUMMER STATS
write.table(summary.df,paste(pathResults,Sys.Date(),"BAGERSummary.txt", sep=""),sep='\t')



