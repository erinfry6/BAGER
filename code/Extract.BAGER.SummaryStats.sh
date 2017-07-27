## this script identifies the best model parameters to include in analyzing the evolution of gene expression
## created by Erin Fry (efry@uchicago.edu)
## User should have R installed
## Non-indented lines should be evaluated for modification speficic to user's purpose

	tissue=$1

###########################################################

export path=/Users/lynchlab/Desktop/ErinFry/workflowr/AGER ##full absolute path to main directory

	export pathData=${path}/data/forBAGER
	export pathScripts=${path}/code
	export pathResults=${path}/data/BAGERresults/${tissue}
	export pathCommands=${pathScripts}/commands
	export pathTemp=${pathResults}/temporary
	export pathAncReconResults=${pathResults}/AncStates
	export pathCommands=${pathScripts}/commands


###########################################################

	## Extract reconstruction information using 'Create.AGER.Summary.File.R'
	
	#R --vanilla <Create.BAGER.Summary.File.R $tissue


