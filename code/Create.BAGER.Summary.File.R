## read in tissue argument

tissue<-"ts"

## set paths to directories, be sure to modify your home directory and the Anccestral Reconstruction directory you are analyzing
path="/Users/lynchlab/Desktop/ErinFry/workflowr/AGER/"
pathData=paste(path,"data/forBAGER/",sep="")
pathResults=paste(path,"data/BAGERresults/",tissue,"/",sep="")
pathAncRecon=paste(pathResults,"AncStates/",sep="")

## set the models tested
models<-c("VarRates", "Lambda")

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

## set function for importing the ancestral expression reconstructions, excluding the header
read.AncRecon=function(file, firstrow, lastrow,header=F, sep='\t'){
  temp<-t(read.tcsv(file, sep='\t', header=F)[,firstrow:lastrow])
  colnames(temp)<-temp[1,]
  temp<-temp[-1,]
  return(temp)
}

#The function DistDiv finds the frequency at which two distributions are different from one another
#the first two arguments are the distributions of interest
#the second is the number of bins you would like to divide the data into, default 10

DistDiv<-function(dist1,dist2,nbin=100) {
  dist2<-as.numeric(dist2)
  dist1<-as.numeric(dist1)
  #first, define the bins each distribution will be broken up into
  minimum=(min(dist1, dist2)) #minimum value of both distributions
  maximum=(max(dist1, dist2)) #maximum value of both distributions
  bins <- seq(minimum, maximum, by =(maximum-minimum)/nbin )  #create nbins from the minimum to maximum observed values
  
  #Create a data frame to contain the number of counts from each distribution in each bin
  #the hist(plot=FALSE) function creates a list containing count information in each bin, speficied above
  counts<-as.data.frame(cbind(hist(dist1, plot=FALSE, breaks=bins)$counts,hist(dist2, plot=FALSE, breaks=bins)$counts))
  colnames(counts)<- c("Dist1Counts", "Dist2Counts") #set the column names
  
  #find the number of overlapping counts across all bins
  ##create new column containing the minimum count of the two distributions
  ##this minimum count is equal to half of the overlap between the two in that bin
  counts$overlap<-apply(counts[,1:2],1,min)  #Take the minimum count for each bin
  
  #multiple the overlap by two to equal the percent overlap between the two distributions
  #then divide by the total number of observations to get the proportion overlap between the two distributions
  return(1-(2*sum(counts$overlap))/sum(counts$Dist1Counts,counts$Dist2Counts))    }


######################################################################

## import the gene information to we can include them in our analysis
genenames<-read.csv(paste(pathData,tissue,"_genesincluded.txt",sep=""),header=T, sep='\t')

## import other infomration of those genes from biomart
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(genenames)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol","ensembl_gene_id","chromosome_name"),values=genes,mart= mart)

## because the first gene is duplicated, duplicate the first gene in G_list to stay consistant
G_list<-rbind(G_list[1,],G_list)

## set the nodes of interest
if (tissue=="br"){
  ancHominini="Node-00003"
  ancHomo="Node-00004"
  
} else if (tissue=="cb"){
  ancHominini="Node-00003"
  ancHomo="Node-00004"
  
} else if (tissue=="lv"){
  ancHominini="Node-00003"
  ancHomo="Node-00004"
  
} else if (tissue=="kd"){
  ancHominini="Node-00003"
  ancHomo="Node-00004"
  
} else if (tissue=="ts"){
  ancHominini="Node-00002"
  ancHomo="Node-00003"
  
} else if (tissue=="ht"){
  ancHominini="Node-00003"
  ancHomo="Node-00004"
  
} else {
  print("this tissue is not known") }

######################################################################
## calculate statistics for the AGERs at the ancHomo and ancHomonini nodes


modelchoice<-t(read.csv(paste(pathResults,"modelchoice.txt", sep=""),header=F, sep="\t"))
colnames(modelchoice)<-modelchoice[1,]
modelchoice<-as.data.frame(modelchoice[-1,])

listcsv<-modelchoice$gene.number

## find which row the iteration information begins on for this tissue's tree
finding.information.about.file<-(read.tcsv(paste(pathAncRecon,models[length(models)],"/",listcsv[1],sep=""), sep='\t'))
it.begin<-which(colnames(finding.information.about.file)=="Itter")

## set the total number of rows you expect to have in each file so the code will warn you if reconstructions failed
expectedrows=ncol(finding.information.about.file)-it.begin


## set empty vectors to contain statistics
foldSD<-vector()
percent.divergent<-vector()
BayesianPostProbofDivergence<-vector()
MedianAncHomo<-vector() # the median reconstructed value of human
MedianAncHominini<-vector() # the median reconstructed value of HC


## for each gene
for (i in 1:nrow(modelchoice)){

  ## load the model choice
choice<-modelchoice$model.choice[i]

## for each gene in the list
  gene<-read.AncRecon(paste(pathAncRecon,choice,"/",listcsv[i],sep=""), firstrow = it.begin,lastrow=(it.begin+expectedrows), sep='\t') # read in the output of BayesTraits
  if (nrow(gene)!=expectedrows){
    warning(paste("The reconstruction of gene number", i, "failed. Check file"))
  }
  
  ## calculate the fold difference increase in the hominini reconstruction
  foldSD[i]<-sd(gene[,which(colnames(gene)==paste(ancHominini, " - 1",sep=""))])/sd(gene[,which(colnames(gene)==paste(ancHomo, " - 1",sep=""))])

  ## caclculate the percent divergence of the two distributions; this is not the preffered method for identifying shifts
  percent.divergent[i]<-DistDiv(gene[,which(colnames(gene)==paste(ancHomo, " - 1",sep=""))], gene[,which(colnames(gene)==paste(ancHominini, " - 1",sep=""))])
  
  ## the preferred method is the Bayesian Posterior Probability of Divergence
  diff<-as.numeric(gene[,which(colnames(gene)==paste(ancHomo, " - 1",sep=""))])-as.numeric(gene[,which(colnames(gene)==paste(ancHominini, " - 1",sep=""))])
  BayesianPostProbofDivergence[i]<-abs(max(1-(length(which(diff>0))/(expectedrows-1)), (length(which(diff>0))/(expectedrows-1))))
  
  ## find the median reconstructed values
  MedianAncHominini[i]<-median(as.numeric(gene[,which(colnames(gene)==paste(ancHominini, " - 1",sep=""))]))
  MedianAncHomo[i]<-median(as.numeric(gene[,which(colnames(gene)==paste(ancHomo, " - 1",sep=""))]))
  
  }

## sometimes there are multiple ensembl id's per gene name, and this will misalign gene labels with the reconstruction information
dup.ensembl<-which(rownames(genenames) %in% G_list$ensembl_gene_id ==FALSE)


## combine gene information, divergence data, convergence data, and means and confidence intervals into one dataframe
Summary<-as.data.frame(cbind(listcsv[-dup.ensembl],modelchoice$model.choice[-dup.ensembl],G_list, BayesianPostProbofDivergence[-dup.ensembl],foldSD[-dup.ensembl],MedianAncHominini[-dup.ensembl],MedianAncHomo[-dup.ensembl],percent.divergent[-dup.ensembl]))
colnames(Summary)<-c("listcsv","modelchoice","hgnc_symbol", "ensembl_gene_id", "chromosome_name", "BayesianPostProbofDivergence", "foldSD", "MedianAncHominini", "MedianAncHomo", "percent.divergent")

## save data
write.table(Summary,paste(pathResults,Sys.Date(),"BAGERSummary.txt", sep=""),sep='\t')

head(Summary)
}


