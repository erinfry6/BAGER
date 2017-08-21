
## Erin Fry
## August 21 2017

## this script consolidates all significant GO enrichments (FDR<0.05) in all lineages after running IdentifyShiftsInAllLineages.Rmd
## User must modify the base path
## The other basic options to change might be the models to evaluate if different than the recommended models

	
###########################################################

	## read in tissue code and start and end gene numbers

	tissue=$1

###########################################################

	## set directory paths

export path=/Users/lynchlab/Desktop/ErinFry/workflowr/AGER ##full absolute path to main directory

	export pathData=${path}/data/forBAGER
	export pathScripts=${path}/code
	export pathResults=${path}/data/BAGERresults/${tissue}
	export pathCommands=${pathScripts}/commands
	export pathEnrichments=${pathResults}/Enrichments
	


###########################################################

## create results file
## first give it a header with column names corresponding to results below for easier manipulation

echo -e 'geneset	description	link	C	O	E	R	PValue	FDR	overlapGene	OverlapGene_UserID' > $pathEnrichments/enrichmentresults.txt

## explanation of what the file contains
echo These are the significant enrichments FDR less than 5 percent for $tissue in all primate lineages >> $pathEnrichments/enrichmentresults.txt
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ >> $pathEnrichments/enrichmentresults.txt

cd $pathEnrichments

for enrich in Project*/enrichment_results*
	do
	echo $enrich has the following enrichments >> $pathEnrichments/enrichmentresults.txt
	cat $enrich >> $pathEnrichments/enrichmentresults.txt

	echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ >> $pathEnrichments/enrichmentresults.txt

done