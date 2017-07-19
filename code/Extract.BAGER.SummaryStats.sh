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

	## for the following code to work, you must have a gene 1. duplicate gene #2
	## define the models to be tested
	models="VarRates Lambda"

	for m in $models
	do

	if [ -e ${pathSSSResults}/$m/gene1.txt ]; then
   	echo gene 1 already duplicated
    else
	cp -r ${pathSSSResults}/$m/gene2.txt ${pathSSSResults}/$m/gene1.txt
    fi

	done

###########################################################

	## Extract reconstruction information using 'Create.AGER.Summary.File.R'
	
	R --vanilla <Create.AGER.Summary.File.R $tissue


