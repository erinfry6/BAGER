
## Erin Fry
## August 22 2017

## this script create demo files to use BAGER from the frontal cortex analysis


######################################################################

## set tissue code

tissue<-"br"

## set paths to directories, be sure to modify your home directory and the Anccestral Reconstruction directory you are analyzing
path="/Users/lynchlab/Desktop/ErinFry/workflowr/AGER/"
pathData=paste(path,"data/forBAGER/",sep="")
pathResults=paste(path,"data/BAGERresults/",tissue,"/",sep="")
pathDemo=paste(path,"demo",sep="")


## read in shift and reconstructed genes
shiftgenes<-read.table(paste(pathResults, "2017-08-22ancHominini_ancHomo90Cutoff.txt", sep=""))
reconsgenes<-read.table(paste(pathResults, "2017-08-22ancHominini_ancHomoReconstructed.txt", sep=""))

## define the genes to keep for demo
Brawandshifts<-shiftgenes[shiftgenes$hgnc_symbol %in% c("LIX1", "CENPT", "SYT15", "THBS4"),]
othershifts<-shiftgenes[1:4,]
nonshifts<-reconsgenes[1:4,] 
## combine into one df
demoresults<-rbind(Brawandshifts, othershifts, nonshifts)

## read in orginal expression data
exp_data<-read.table(paste(pathData, "br_genesincluded.txt", sep=""))
colnames(exp_data)<-exp_data[1,]
exp_data<-exp_data[-1,]

## take expression data for demo genes
demo_exp_data<-exp_data[exp_data$EnsemblID %in% demoresults$ensembl_gene_id,]

## save files
write.table(demo_exp_data, paste(pathDemo, "/br_genesincluded.txt", sep=""), sep="\t", row.names = FALSE)
write.table(t(demo_exp_data)[-1,], paste(pathDemo, "/br_exp.txt", sep=""), sep="\t", col.names = FALSE)

