# Ancestral Transcriptome Reconstruction

The following scripts can be used to reconstruct ancestral posterior probability distributions of gene expression at internal nodes within a phylogenetic tree. These scripts and pipeline were written by Erin Fry (efry@uchicago.edu) in the Lynch Laboratory, Department of Human Genetics, University of Chicago.

They will:

1) Reconstruct ancestral gene expression levels using BayesTrait’s continuous trait, random walk MCMC algorithm, transcriptome data from extant species, and the known species phylogeny with distances.

2) Identify genes with expression shifts in the lineage of interest by comparing the posterior probability distributions of the ancestral reconstructions

## Set up directories

Before beginning, create a home directory for the pipeline that contains the following subdirectories

				home/data/
					
							forBAGER
							
							BAGERresults
				
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

 -  Please see the instructions at the top of each individual script for modification requirements.

 - If you would like to modify the models ran in BayesTraits from the default models in this pipeline, you will need to specific which models in`create_model_files.sh` , `identify_best_model.sh` , `ID_best_model.R`, and `Create.BAGER.Summary.File.R`.

 
## Run the Bayesian Ancestral Transcriptome Reconstruction Scripts


#### 1) create_model_files.sh - Run the MCMC chain under desired models. The default recommended models are coded in, but if you would like to change or add models, you will have to do so.
I typically run the code one tissue at a time, broken into 4 or 5 chunks specified by begingene and endgene to decrease run time per tissue.

```
./create_model_files.sh tissuecode begingene endgene
```


#### 2) ID_best_model.R - Identify the model with the largest log marginal likelihood. Saves model choice in the results as modelchoice.txt in the BAGERresults/$tissuecode directory

```
R --vanilla <ID_best_model.R
```


#### 3) Create.BAGER.Summary.File.R - Collect Summary statistics for the BAGERs into one Summary file.

```
.R --vanilla <Create.BAGER.Summary.File.R tissuecode
```

#### 4) BAGERAnalysis.Rmd - Analyzes Ancestral Transcriptome Reconstructions to identify genes with expression shifts. Best used in R studio.



## Simulate gene expression evolution across your tree.

**Using the `SimulateAGER.Rmd` file, you may simulate data that matches your data's expression patterns.**


## Troubleshooting on test example files

To run these scripts on the example test files, copy the files in scripts/test to the data directory and follow the above steps.

### written by Erin Fry
### Last modified: April 24 2017