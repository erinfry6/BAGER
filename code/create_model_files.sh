## this script runs the BayesTraits program given a set of commands for the genes indicated in the for loop
## created by Erin Fry (efry@uchicago.edu)
## Non-indented lines should be evaluated for modification specific to the user's purpose

###########################################################
	## read in which tissue

	tissue=$1

###########################################################

	## set directory paths

export path=/Users/lynchlab/Desktop/ErinFry/workflowr/AGER ##full absolute path to main directory

	export pathData=${path}/data/forBADGER
	export pathScripts=${path}/code
	export pathResults=${path}/data/BADGERresults/${tissue}
	export pathTemp=${pathResults}/temporary
	export pathModelResults=${pathResults}/Models
	export pathSSSResults=${pathResults}/SSS
	export pathCommands=${pathScripts}/commands
	
	## define pathways of files to be created and utilized around during analysis
	export scriptversion=$tissue  ## be sure to save files for their tissue as to not confuse
	export singleexpression=${pathTemp}/singlegene_$scriptversion.txt
	export model=${pathTemp}/Model_$scriptversion.bin
	export commandfile=${pathCommands}/step1command_$scriptversion.txt
	
###########################################################
	
	## make directories if they do not already exist
	
	if [ -e ${pathResults} ]; then
   	echo 'Results dir already here'
    else
    mkdir ${pathResults}
    fi 
	
	if [ -e ${pathTemp} ]; then
   	echo 'Temporary dir already here'
    else
    mkdir ${pathTemp}
    fi 
	
	if [ -e ${pathCommands} ]; then
   	echo 'Command dir already here'
    else
    mkdir ${pathCommands}
    fi
    
	## make directory to store the results from the chain
	if [ -e ${pathMCMCResults} ]; then
   	echo 'MCMC Results dir already here'
    else
    mkdir ${pathMCMCResults}
    fi
    
	## make directory to store the model from the chain
	if [ -e ${pathModelResults} ]; then
   	echo 'MCMC Results dir already here'
    else
    mkdir ${pathModelResults} ${pathModelResults}/kappa ${pathModelResults}/delta ${pathModelResults}/none ${pathModelResults}/kd
    fi
    
    ## make directory to store the results from the stepping stone sampler
	if [ -e ${pathSSSResults} ]; then
   	echo 'SSS Results dir already here'
    else
    mkdir ${pathSSSResults} ${pathSSSResults}/kappa ${pathSSSResults}/delta ${pathSSSResults}/none ${pathSSSResults}/kd
    fi

###########################################################

	## Read in files for tissue of interest
	Expressiondata=${pathData}/${tissue}_exp.txt
	
	## number of genes evaluating
	#NumGenes=$(awk '{print NF}' $Expressiondata | tail -n 1)
	NumGenes=3
	
	echo Evaluating ${NumGenes} genes in ${tissue}

	## define tree file
	tree=${pathData}/${tissue}_tree.tree

###########################################################

	## Run MCMC under each model, saving the summary and log likelihood outputs

	## The descriptions of each section are explained in the following lines with labels to the right of the commands in the scripts, but generally:
	## for each model to be tested, make the appropriate command file, then for each of the genes in your expression data file,
	## create two files summarizing the MCMC run for that trait's evolution and find the likelihood of the model given the data
	
	## models to test: name each of the models you are testing
	## be sure that a directory was created for each of the models (of the same name) in ${pathModelResults}
	
	## general commands: you may modify these general commands referring to the properties of the MCMC, see the BayesTraits manual for more details
	## It is not recommended to remove or modify the general commands that are written here. '4' and '2' indicated a continuous trait random walk model.
	## using the stepping stone sampler is strongly recommended (stones command)
	## Any additional commands you add to each of the models tested should be added here
	
	## model-specific commands: these commands are the lines necessary to run each of the models you wish to test
	
	## second general commands: these commands are necessary to save model files and run the program, do not modify
	
	## create models for each gene: for each gene in 2 through the number of genes in the expression data file, run BayesTraits and save the necessary files
	## each gene will be modeled under each model specified in models to test, and saved accordingly

for c in kd kappa delta none; 											## models to test
	do

echo '4
2
Iterations 1010000
Burnin 10000
stones 100 10000' > ${commandfile}  									## general commands

	if [[ $c == none ]]; then
   	echo creating 'none' model files for ${tissue} 						## model-specific commands

	elif [[ $c == kappa ]]; then
	echo creating kappa model files for ${tissue}
	echo kappa >> ${commandfile}

	elif [[ $c == kd ]]; then
	echo creating kappa and delta model files for ${tissue}
	echo 'kappa
	delta' >> ${commandfile}
	
	else 
	echo creating delta model files for ${tissue}
	echo delta >> ${commandfile}
	fi

	echo SaveModels $model >> ${commandfile} 							## second general commands
	echo run >> ${commandfile}

		

	for ((a=2; a<=$NumGenes; a++))										## create models for each gene
		do
		awk -v a="$a" '{print $1,$a}' ${Expressiondata} > ${singleexpression}

		./../BayesTraits/BayesTraitsV3 ${tree} ${singleexpression} <${commandfile} > ${pathTemp}/MCMC.txt
		cp ${singleexpression}.Stones.txt ${pathSSSResults}/${c}/gene$a.txt
		cp ${model} ${pathModelResults}/${c}/gene$a.bin

		done		


done
		