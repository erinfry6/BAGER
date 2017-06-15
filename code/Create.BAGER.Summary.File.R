## read in tissue argument

tissue<-commandArgs()

## set paths to directories, be sure to modify your home directory and the Anccestral Reconstruction directory you are analyzing
path="/Users/lynchlab/Desktop/ErinFry/workflowr/AGER/"
pathData=paste(path,"data/forBAGER/",sep="")
pathResults=paste(path,"data/BAGERresults/",tissue,"/",sep="")
pathAncRecon=paste(pathResults,"AncRecon/",sep="")

library(dplyr)

######################################################################

setwd(pathAncRecon)
options(stringsAsFactors = FALSE)
ldf <- list() # creates a list
listcsv<-paste("gene",as.character((order(dir(pattern = "*.txt")))),".txt", sep="")

## for each gene, read into ldf collection of dataframes
## but first set the total number of rows you expect to have in each file so the code will warn you if reconstructions failed
expectedrows<-1001

for (k in 1:length(listcsv)){ 
  ldf[[k]]<- read.csv(listcsv[k], sep='\t') # read files in listcsv into the ldf list
  if (nrow(ldf[[k]])!=expectedrows){
    warning(paste("The reconstruction of gene number", k, "failed. Check file"))
  }
}

## import the gene information to we can include them in our analysis
setwd(path)
genenames<-read.csv("GeneNamesandcodes.txt",header=T, sep='\t')

######################################################################


### Write dist div function

#The function DistDiv finds the frequency at which two distributions are different from one another
#'dylpr' is a required package for this function
#the first two arguments are the distributions of interest
#the second is the number of bins you would like to divide the data into, default 100

DistDiv<-function(dist1,dist2,nbin=100) {
  
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


## Summarize the AGERs
## Find AGERs with succesful reconstructions
## Genes that succesfully reconstructed will have similar posterior variations of the two ancestral node reconstructions. Thus, to eliminate successfully reconstructed genes, you calculate the fold difference in standard deivation between the two reconstructions.


## the HC ancestral reconstruction generally has a higher standard deviation. If it is too much larger than the ancHuman the chain failed to reconstruct

## find the fold difference in standard deviation between the two reconstructions

foldSD<-vector() ## empty vector to store results
for (i in 1:length(ldf)){
  foldSD[i]<-sd(ldf[[i]]$Est.AncHominini...1[1:1000])/sd(ldf[[i]]$Est.AncHomo...1[1:1000])
}


### Calculate divergence between the two reconstructions using two methods- Percent Divergence between the two Post Prob Distn's of Expression and by the Bayesian Posterior Probability of Divergence
  
percent.divergent<-vector(length=length(ldf)) #create empty vectors to contain percent divergence of each gene

for (i in 1:length(ldf)){  #for each of these genes
  gene<- ldf[[i]][-nrow(ldf[[i]]),]  #create a new file for only that gene
  
  #run the function distdiv for the human and H-C post prob distributions
  percent.divergent[i]<-DistDiv(gene$Est.AncHomo...1, gene$Est.AncHominini...1)
}


## Calculate the Bayesian Posterior Probability of Divergence (define in introduction)
BayesianPostProbofDivergence<-vector(length=length(ldf))

for (i in 1:length(ldf)){  #for each of these genes
  gene<- ldf[[i]][-nrow(ldf[[i]]),]  #create a new file for only that gene
  diff<-gene$Est.AncHomo...1-gene$Est.AncHominini...1
  
  #find the proportion that are greater than 0
  BayesianPostProbofDivergence[i]<-abs(max(1-(length(which(diff>0))/(expectedrows-1)), (length(which(diff>0))/(expectedrows-1))))
}



### Calculate the Posterior Median of the distributions

#find the mean of the ancestral estimations, and for fun, the lower and upper confidence intervals
MedianAncHomo<-vector(length=length(ldf)) # the median of human
MedianAncHominini<-vector(length=length(ldf)) # the median of HC

for (i in 1:length(ldf)){
  MedianAncHominini[i]<-median(ldf[[i]][-nrow(ldf[[i]]),]$Est.AncHominini...1)
  
  MedianAncHomo[i]<-median(ldf[[i]][-nrow(ldf[[i]]),]$Est.AncHomo...1 )
  
}


######################################################################

### Combine into Summary file

## combine gene information, divergence data, convergence data, and means and confidence intervals into one dataframe
Summary<-as.data.frame(cbind(genenames[1:length(ldf),], percent.divergent,BayesianPostProbofDivergence,foldSD,MedianAncHominini,MedianAncHomo))

setwd(pathResults)

## save data
write.table(Summary,"AGERSummary.txt",sep='\t')

head(Summary)



