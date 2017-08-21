## this script reconstructs the ancestral state, or transcription level, of your transcriptome
## created by Erin Fry (efry@uchicago.edu)
## Non-indented lines should be evaluated for modification specific to the user's purpose

###########################################################
	## read in which tissue

	tissue=$1

###########################################################

	## set directory paths

export path=/Users/lynchlab/Desktop/ErinFry/workflowr/AGER ##full absolute path to main directory

	export pathData=${path}/data/forBAGER
	export pathScripts=${path}/code
	export pathResults=${path}/data/BAGERresults/${tissue}
	export pathTemp=${pathResults}/temporary
	export pathModelResults=${pathResults}/Models
	export pathSSSResults=${pathResults}/ARSSS
	export pathCommands=${pathScripts}/commands
	export pathRecon=${pathResults}/AncRecon
	
	## define pathways of files to be created and utilized around during analysis
	export scriptversion=$tissue  ## be sure to save files for their tissue as to not confuse
	export singleexpression=${pathTemp}/singlegene_$scriptversion.txt
	export model=${pathTemp}/Model_$scriptversion.bin
	export commandfile=${pathCommands}/step2command_$scriptversion.txt

###########################################################

	## make directory to store the reconstructions
	
	if [ -e ${pathRecon} ]; then
   	echo 'AncRecon dir already here'
    else
    mkdir ${pathRecon}
    fi
 	
 	## make directory to store the results from the stepping stone sampler
	if [ -e ${pathSSSResults} ]; then
   	echo 'SSS Results dir already here'
    else
    mkdir ${pathSSSResults}
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

	## Reconstruct the ancestral gene expression of each gene in the dataset, saving the likelhhood and Bayesian Posterior Distributions of the Ancestral States

	## The descriptions of each section are explained in the following lines with labels to the right of the commands in the scripts, but generally:
	## for each gene, make the appropriate command file based on the model that best fits the gene's expression
	## run the MCMC and save the reconstructed inferred expression values
	
	## genes to test: go through each gene
	
	## general commands and ancestral nodes: Any additional commands you add to each of the models tested should be added here
	## specifies which ancestors to reconstruct: instructions in BayesTraits V3 manual
	## You may modify these general commands referring to the properties of the MCMC, see the BayesTraits manual for more details
	## It is not recommended to remove or modify most of the general commands that are written here. '4' and '2' indicated a continuous trait random walk model.
	## Using the stepping stone sampler is strongly recommended (stones command)
	
	## model-specific commands and reconstruction: these commands are the lines necessary to run each of the models you wish to test, modify them using echo >> ${commandfile}
	## last line in each must end in 'run'
	## each gene will be reconstructed under its best fit model, and saved



###########################################################

## for loop goes through each of the genes specified in { .. }
## choice of 1=Kappa, 2= KD, 3=Delta, 4=none
	## load in the model choice
	## make a temporary file to contain only gene expression from the one gene
	## copy the model file created in the first step and run BayesTraits again, informed by the model file to reconstruct ancestral transcriptional state
	## copy the stepping stone output to save likelihood information about the chain

for ((a=2; a<=$NumGenes; a++))													## genes to test
	do
	
echo '7
2
VarRates
Iterations 110000000
Burnin 10000000
stones 100 10000
AddTag Tag-PointA hsa_br_M_2 hsa_br_M_3 hsa_br_F_1 
AddTag Tag-PointB hsa_br_M_2 hsa_br_M_3 hsa_br_F_1 ptr_br_M_3 ptr_br_M_2 ptr_br_M_5 ptr_br_M_1 ptr_br_M_4 ptr_br_F_1 ppa_br_M_1 ppa_br_F_1 ppa_br_F_2 
AddNode AncHomo Tag-PointA
AddNode AncHominini Tag-PointB' >> ${commandfile}
echo LoadModels ${pathTemp}/model$scriptversion.bin >> ${commandfile} 			## general commands and ancestral nodes
	
	choice=$(awk -v a="$a" '{print $a}' ${pathResults}/modelchoice.txt) || exit 1
	
	awk -v a="$a" '{print $1,$a}' ${pathData}/${Expressiondata} > ${singleexpression} || exit 1
	
if [[ $choice == 4 ]]; then

echo run >> ${commandfile}														## model-specific commands and reconstruction
	
cp ${pathModelResults}/none/gene$a.bin ${pathTemp}/model$scriptversion.bin
./../BayesTraits/BayesTraitsV3 ${tree} ${singleexpression} <${commandfile} | awk 'NR >=91' > ${pathRecon}/gene$a.txt || exit 1
	

elif [[ $choice == 3 ]]; then

echo 'delta'																	## model-specific commands and reconstruction
run' >> ${commandfile}
	
cp ${pathModelResults}/delta/gene$a.bin ${pathTemp}/model$scriptversion.bin
./../BayesTraits/BayesTraitsV3 ${tree} ${singleexpression} <${commandfile} | awk 'NR >=92' > ${pathRecon}/gene$a.txt || exit 1

elif [[ $choice == 2 ]]; then

echo 'kappa
delta
run' >> ${commandfile}															## model-specific commands and reconstruction
	
cp ${pathModelResults}/kd/gene$a.bin ${pathTemp}/model$scriptversion.bin
./../BayesTraits/BayesTraitsV3 ${tree} ${singleexpression} <${commandfile} | awk 'NR >=93' > ${pathRecon}/gene$a.txt || exit 1
	
else 

echo 'kappa
run' >> ${commandfile}															## model-specific commands and reconstruction
	
cp ${pathModelResults}/kappa/gene$a.bin ${pathTemp}/model$scriptversion.bin
./../BayesTraits/BayesTraitsV3 ${tree} ${singleexpression} <${commandfile} | awk 'NR >=92' > ${pathRecon}/gene$a.txt || exit 1
	
	fi
	
	cp ${singleexpression}.Stones.txt ${pathSSSResults}/gene$a.txt || exit 1

done
