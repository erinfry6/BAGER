# Ancestral Transcriptome Reconstruction

The following scripts can be used to reconstruct ancestral posterior probability distributions of gene expression at internal nodes within a phylogenetic tree. These scripts and pipeline were written by Erin Fry (efry@uchicago.edu) in the Lynch Laboratory, Department of Human Genetics, University of Chicago.

They will:

1) Reconstruct ancestral gene expression levels using BayesTrait’s continuous trait, random walk MCMC algorithm, transcriptome data from extant species, and the known species phylogeny with distances.

2) Identify genes with expression shifts in the lineage of interest by comparing the posterior probability distributions of the ancestral reconstructions

## Set up directories

Before beginning, create a home directory for the pipeline that contains the following subdirectories

				home/data/
					
							forBADGER
							
							BADGERresults
				
				home/code
				
				home/BayesTraits
				
				
Place the contents of this repository in the code folder, except for the BayesTraitsV3. Place this in the home/BayesTraits folder.

[BayesTraits Manual](http://www.evolution.rdg.ac.uk/BayesTraitsV3/Files/BayesTraitsV3.Manual.pdf)


## Input files Formats

The following input files are required for this analysis:

 - tissuecode_exp.txt: a tab delimited file containing expression data formatted according to the BayesTraits Manual: Expression data file with no column names. First column: names of the samples that coordinate with the Nexus tree file. All subsequent columns are the expression data for each gene. Be sure to keep another document that notes the gene in each column, as well as gene numbers, beginning at 2.

 - tissuecode_tree.tree: Tree.tree is an ultrametric nexus formatted phylogeny formatted according to the BayesTraits Manual.

* tissuecode is the name of the tissue you are analyzing that will be an input for the scripts below *


## Modify the scripts

 - The top of each script (all .sh and .R files) must be modified to contain the proper home directory path

 - The gene numbering will begin at 2, not 1. 

 - If you would like to modify the models ran in BayesTraits, you will need to manually change the command sections of `create_model_files.sh` and `ancestral_reconstruction.sh`. 

 - To specify which ancestral nodes are reconstructed, modify `ancestral_reconstruction.sh`'s command section, which is at the bottom in the for loop. Follow the instructions in AGERinstructions.txt

## Run the Bayesian Ancestral Transcriptome Reconstruction Scripts


#### 1) create_model_files.sh - Run the MCMC chain under 4 models (without tree transformation parameters, with kappa, with delta, and with both kappa and delta).
In our experience, incorporating both Kappa and Delta (see BayesTraits manual) decreases the variance of reconstructions and increases the log likelihood of the chain.

```
./create_model_files.sh tissuecode
```


#### 2) identify_best_model.sh - Identify the model with the largest log marginal likelihood. Saves model choice in the results as modelchoice.txt

```
./identify_best_model.sh 
```


#### 3) ancestral_reconstruction.sh - Reconstructs the posterior probability distribution of ancestral states at specified internal nodes using the best model for the gene's evolution. 
   These nodes must be specified in the CREATING COMMAND FILES section. See BayesTraits Manual for instructions.

```
./ancestral_reconstruction.sh Expression_Data.txt SampleTree.tree 4 2
```
   
_All four inputs following the bash script should be the exact same as in step 1._


#### 4) Extract.AGER.SummaryStats.sh - Collect Summary statistics for the AGERs into one Summary file using Create.AGER.Summary.File.R .

```
./Extract.BADGER.SummaryStat.sh
```

#### 5) AGERAnalysis.Rmd - Analyzes Ancestral Transcriptome Reconstructions to identify genes with expression shifts. Best used in R studio.



## Simulate gene expression evolution across your tree.

**Using the `SimulateAGER.Rmd` file, you may simulate data that matches your data's expression patterns.**


## Troubleshooting on test example files

To run these scripts on the example test files, copy the files in scripts/test to the data directory and follow the above steps.

### written by Erin Fry
### Last modified: April 24 2017