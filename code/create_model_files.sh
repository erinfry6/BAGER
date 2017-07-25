## this script runs the BayesTraits program given a set of commands for the genes indicated in the for loop
## created by Erin Fry (efry@uchicago.edu)
## Non-indented lines should be evaluated for modification specific to the user's purpose

###########################################################
	## read in which tissue

	tissue=$1
	begingene=$2
	endgene=$3

###########################################################

	## set directory paths

export path=/Users/lynchlab/Desktop/ErinFry/workflowr/AGER ##full absolute path to main directory

	export pathData=${path}/data/forBAGER
	export pathScripts=${path}/code
	export pathResults=${path}/data/BAGERresults/${tissue}
	export pathCommands=${pathScripts}/commands
	export pathTemp=${pathResults}/temporary

	## make directories if they do not already exist
	
	if [ -e ${pathResults} ]; then
   	echo 'Results dir already here'
    else
    mkdir ${pathResults}
    fi 
	
	
	if [ -e ${pathCommands} ]; then
   	echo 'Command dir already here'
    else
    mkdir ${pathCommands}
    fi	
    
	if [ -e ${pathTemp} ]; then
   	echo 'Command dir already here'
    else
    mkdir ${pathTemp}
    fi	

###########################################################

	
	## define the models to be tested
	models="VarRates"
	
	## make a list of .txt files to be created from the model, these should blank in the *.txt______.txt
	## if using VarRates with stepping stone sampler these will be the files created
	filestomake="AncStates Log Schedule Stones VarRates"
	
	for f in $filestomake 
	do
		mkdir ${pathResults}/$f
		for m in $models
		do
			mkdir ${pathResults}/$f/$m
		done
	done

	## make directories for the output tree files for each model
	mkdir ${pathResults}/Output_trees
	for m in $models
		do
			mkdir ${pathResults}/Output_trees/$m
		done
	
	## make directories for the output MCMC files for each model
	mkdir ${pathResults}/MCMC
	for m in $models
		do
			mkdir ${pathResults}/MCMC/$m
		done


###########################################################
	
	## define pathways of files to be created and utilized during analysis
	export scriptversion=$begingene  ## be sure to save files for their tissue as to not confuse, or this can be used to run multiple create_model_files at once
	export singleexpression=${pathResults}/temporary/singlegene_$scriptversion.txt
	export commandfile=${pathCommands}/step1command_$scriptversion.txt
	
###########################################################

	## Read in files for tissue of interest
	Expressiondata=${pathData}/${tissue}_exp.txt
	
	## number of genes evaluating
	NumGenes=$(awk '{print NF}' $Expressiondata | tail -n 1)
	#NumGenes=10
	echo Evaluating ${NumGenes} genes in ${tissue}

	## define tree input file
	tree=${pathData}/${tissue}_tree.tree

###########################################################

	## Run MCMC under each model, saving all output files

	## The descriptions of each section are explained in the following lines with labels to the right of the commands in the scripts, but generally:
	## for each model to be tested, make the appropriate command file, then for each of the genes in your expression data file,
	## run BayesTraits to infer evolutionary rate parameters of the gene's expression across the tree, including ancestral states and rates along branches
	
	## general commands: you may modify these general commands referring to the properties of the MCMC, see the BayesTraits manual for more details
	## 12 is **** special****
	## VarRates tests the variable rates model, a model that breaks the Brownian Motion assumption
	## Using the stepping stone sampler is strongly recommended (stones command)
	## Any additional commands you add to each of the models tested should be added here
	
	## model-specific commands: these commands are the lines necessary to run each of the models you wish to test
	
	## run the program, do not modify
	
	## run model for each gene: for each gene in 2 through the number of genes in the expression data file (unless otherwise noted above), run BayesTraits and save output files
	## each gene will be modeled under each model specified in models to test, and saved accordingly

for m in $models 
do

echo '12
VarRates
stones 100 10000' > ${commandfile}  								## general commands

if [[ $m == VarRates ]]; then										## model-specific commands
echo creating variable rates model files for ${tissue} 
	
elif [[ $m == Lambda ]]; then
echo creating variable rates and lambda model files for ${tissue} 
echo lambda >> ${commandfile}


else 
echo bad model call
fi


echo run >> ${commandfile} 										## second general commands

		

	for ((a=$begingene; a<=$endgene; a++))										## run model for each gene
		do
		
		if [ -e ${pathResults}/AncStates/$m/gene$a.txt ]; then
   		echo already created $m model file for gene $a
    	else
    	
		awk -v a="$a" '{print $1,$a}' ${Expressiondata} > ${singleexpression}

		./../BayesTraits/BP3.1 ${tree} ${singleexpression} <${commandfile} > ${pathTemp}/MCMC$scriptversion.txt
		
		
		fi
		
for f in $filestomake 
do
cp ${singleexpression}.$f.txt ${pathResults}/$f/$m/gene$a.txt
done
		
cp ${pathTemp}/MCMC$scriptversion.txt ${pathResults}/MCMC/$m/gene$a.txt
cp ${singleexpression}.Output.trees ${pathResults}/Output_trees/$m/gene$a.trees 
	

done		


done
		