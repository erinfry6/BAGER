

## Erin Fry
## August 1 2017

## this script finds GO enrichments in all lineages after running IdentifyShiftsInAllLineages.Rmd
## Non-indented lines should be evaluated for modification specific to the user's purpose
## User MUST modify the base path
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

echo These are the enrichments for $tissue in all lineages > $pathEnrichments/enrichmentresults.txt
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ >> $pathEnrichments/enrichmentresults.txt

cd $pathEnrichments

for enrich in Project*/enrichment_results*
do
echo $enrich has the following enrichments >> $pathEnrichments/enrichmentresults.txt
cat $enrich >> $pathEnrichments/enrichmentresults.txt

echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ >> $pathEnrichments/enrichmentresults.txt

done