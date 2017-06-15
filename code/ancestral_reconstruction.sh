## this script reconstructs the ancestral state, or transcription level, of your transcriptome
## created by Erin Fry (efry@uchicago.edu)
## Non-indented lines should be evaluated for modification specific to user's purpose

###########################################################

	## set directory paths

export path=/Users/lynchlab/Desktop/ErinFry/ReconAncNeoTranscriptomes/BrainConstitiutive/BTReconstruct ##full absolute path to main directory

	export pathData=${path}/data
	export pathScripts=${path}/scripts
	export pathResults=${path}/resultsFour
	export pathTemp=${pathResults}/temporary
	export pathModelResults=${pathResults}/Models
	export pathCommands=${pathScripts}/commands
	export pathSSSResults=${pathResults}/ARSSS
	export pathRecon=${pathResults}/AncRecon

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
	## the first command is the file with all gene expression data, the second is the phylogenetic tree with distances
	## the next two are the two commands needed to run BayesTraits (for continuous walk MCMC, 4 and 2, respectively)

	Expressiondata=$1
	tree=$2
	command1=$3
	command2=$4

## if running multiple scripts at once, make sure to modify scriptversion to be a different number in each script to avoid creating the same numbered files

export scriptversion=1  ## modify this if running multiple files at once
	export singleexpression=${pathTemp}/singlegene$scriptversion.txt
	export MCMC=${pathTemp}/MCMC$scriptversion.txt
	export commandfiletemplate=${pathCommands}/step2commandtemplate$scriptversion.txt
	export commandfile=${pathCommands}/step2command$scriptversion.txt
	
###########################################################

## CREATING COMMAND FILES

	## creates the command files to use or not use each evolutionary rate parameter

	echo $command1 > ${commandfiletemplate}
	echo $command2 >> ${commandfiletemplate}

## specifies how many iterations, the burnin period, and the number of stones to sample
## it is important to read AGERinstructions.txt and the BayesTraits manual and modify the label and node numbers to specify which ancestral node to reconstruct

echo 'Iterations 1010000
AddTag Tag-PointA hsa_br_M_1 hsa_br_M_2 hsa_br_M_3 hsa_br_M_4 hsa_br_F_1 
AddTag Tag-PointB hsa_br_M_1 hsa_br_M_2 hsa_br_M_3 hsa_br_M_4 hsa_br_F_1 ptr_br_M_3 ptr_br_M_2 ptr_br_M_5 ptr_br_M_1 ptr_br_M_4 ptr_br_F_1 ppa_br_M_1 ppa_br_F_1 ppa_br_F_2 
AddMRCA AncHomo Tag-PointA
AddMRCA AncHominini Tag-PointB
Prior AncState-1 uniform 0 15731
Burnin 10000
stones 100 10000' >> ${commandfiletemplate}
echo LoadModels ${pathTemp}/model$scriptversion.bin >> ${commandfiletemplate}


###########################################################

## for loop goes through each of the genes specified in { .. }
## choice of 1=Kappa, 2= KD, 3=Delta, 4=none
	## load in the model choice
	## make a temporary file to contain only gene expression from the one gene
	## copy the model file created in the first step and run BayesTraits again, informed by the model file to reconstruct ancestral transcriptional state
	## copy the stepping stone output to save likelihood information about the chain

for a in {2..3300}
	do
	
	choice=$(awk -v a="$a" '{print $a}' ${pathResults}/modelchoice.txt)
	
	awk -v a="$a" '{print $1,$a}' ${pathData}/${Expressiondata} > ${singleexpression}
	
	cp ${commandfiletemplate} ${commandfile}
	
	if [[ $choice == 4 ]]; then
		echo run >> ${commandfile}
	
		cp ${pathModelResults}/none/gene$a.bin ${pathTemp}/model$scriptversion.bin
		./../BayesTraitsV2/BayesTraitsV2 ${pathData}/${tree} ${singleexpression} <${commandfile} | awk 'NR >=91' > ${pathRecon}/gene$a.txt
	

	elif [[ $choice == 3 ]]; then
		echo 'delta
		run' >> ${commandfile}
	
		cp ${pathModelResults}/delta/gene$a.bin ${pathTemp}/model$scriptversion.bin
		./../BayesTraitsV2/BayesTraitsV2 ${pathData}/${tree} ${singleexpression} <${commandfile} | awk 'NR >=92' > ${pathRecon}/gene$a.txt
	
	elif [[ $choice == 2 ]]; then
		echo 'kappa
		delta
		run' >> ${commandfile}
	
		cp ${pathModelResults}/kd/gene$a.bin ${pathTemp}/model$scriptversion.bin
		./../BayesTraitsV2/BayesTraitsV2 ${pathData}/${tree} ${singleexpression} <${commandfile} | awk 'NR >=93' > ${pathRecon}/gene$a.txt
	
	else 
		echo 'kappa
		run' >> ${commandfile}
	
		cp ${pathModelResults}/kappa/gene$a.bin ${pathTemp}/model$scriptversion.bin
		./../BayesTraitsV2/BayesTraitsV2 ${pathData}/${tree} ${singleexpression} <${commandfile} | awk 'NR >=92' > ${pathRecon}/gene$a.txt
	fi
	cp ${singleexpression}.log.txt.Stones.txt ${pathSSSResults}/gene$a.txt

done
