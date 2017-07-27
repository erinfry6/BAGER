## this script moves genes up. need to redo the first 17 genes of testis and move the rest

## only tested on ts. 

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

################

	
	## define the models to be tested
	models="VarRates"
	
	## make a list of .txt files to be created from the model, these should blank in the *.txt______.txt
	## if using VarRates with stepping stone sampler these will be the files created
	filestomake="AncStates Log MCMC Schedule Stones VarRates"
	

	
###########################################################

	
	## number of genes evaluating
	echo Moving genes $begingene thru $endgene in ${tissue}


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
		

for ((a=$begingene; a<=$endgene; a++))										## run model for each gene
do
    	
b=$((a - 1))

mv ${pathResults}/Output_trees/$m/gene$a.trees ${pathResults}/Output_trees/$m/gene$b.trees
		
for f in $filestomake 
do
mv ${pathResults}/$f/$m/gene$a.txt ${pathResults}/$f/$m/gene$b.txt
done

done		


done