## this script identifies the best model parameters to include in analyzing the evolution of gene expression
## created by Erin Fry (efry@uchicago.edu)
## User should have R installed
## Non-indented lines should be evaluated for modification speficic to user's purpose

###########################################################

export path=/Users/lynchlab/Desktop/ErinFry/workflowr/AGER ##full absolute path to main directory
	
	export pathData=${path}/data/forBAGER
	export pathScripts=${path}/code
	export pathResults=${path}/data/BAGERresults/${tissue}
	export pathTemporary=${pathResults}/temporary
	export pathCommands=${pathScripts}/commands

###########################################################

	## for the following code to work, you must have a gene 1. for me this is a duplicate of gene #2

	if [ -e ${pathResults}/AncRecon/gene1.txt ]; then
   	echo 'already here'
    else
    cp -r ${pathResults}/AncRecon/gene2.txt ${pathResults}/AncRecon/gene1.txt
    fi


###########################################################

	## Extract reconstruction information using 'Create.AGER.Summary.File.R'
	
	R --vanilla <Create.AGER.Summary.File.R $tissue


