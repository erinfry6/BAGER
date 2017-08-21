# Ancestral Transcriptome Reconstruction

The following scripts can be used to identify gene expression shifts at internal branches of a phylogenetic tree. These scripts and pipeline were written by Erin Fry (efry@uchicago.edu) in the Lynch Laboratory, Department of Human Genetics, University of Chicago.

They will:

1) Reconstruct ancestral gene expression levels using a customized BayesTraits variable rate MCMC algorithm for continuous traits, transcriptome data from extant species, and the known species phylogeny with distances (in time / millions of years).

2) Identify genes with expression shifts across the phylogeny by comparing the posterior probability distributions of the ancestral reconstructions.

## Set up directories

Before beginning, create a home directory for the pipeline that contains the following subdirectories

				home/data/
					
					forBAGER
							
					BAGERresults
				
				home/code
				
				home/BayesTraits
				
				
Place the contents of this repository in the code folder, except for the BP3.1. Place this in the home/BayesTraits folder.

[BayesTraits Manual](http://www.evolution.rdg.ac.uk/BayesTraitsV3/Files/BayesTraitsV3.Manual.pdf)


## Input files Formats

The following input files are required for this analysis and should be places in `home/data/forBAGER`:

 - tissuecode_exp.txt: a tab delimited file containing expression data formatted according to the BayesTraits Manual: Expression data file with no column names. First column: names of the samples that coordinate with the Nexus tree file. All subsequent columns are the expression data for each gene. Be sure to keep another document that notes the gene in each column.

 - tissuecode_tree.tree: tissuecode_tree.tree is an ultrametric nexus formatted phylogeny formatted according to the BayesTraits Manual.
 
 - tissuecode_genesincluded.txt: a tab delimited file containing expression data identical to that in tissuecode_exp.txt, but formatted differently. This file should have column and row names. Columns are the samples and rows are each gene, labeled with its EnsemblID.

* tissuecode is the name of the tissue you are analyzing that will be an input for the scripts below *


## Modify the scripts

 - **Please see the instructions at the top of each individual script for modification requirements.**

 - If you would like to modify the models testeed by BayesTraits from the default models in this pipeline, you will need to specify which models in`create_model_files.sh` , `identify_best_model.sh` , `ID_best_model.R`, and `Create.BAGER.Summary.File.R`.

 
## Run the Bayesian Ancestral Transcriptome Reconstruction Scripts


#### 1) create_model_files.sh - Reconstruct transcriptomes under each model. The default recommended models are coded in, but if you would like to change or add models, you will have to do so.
I typically run the code one tissue at a time, broken into 4 or 5 chunks specified by begingene and endgene numbers to decrease run time per tissue. Begingene and endgene are integers of the beginning and start gene number you would like to reconstruct.

```
./create_model_files.sh tissuecode begingene endgene
```


#### 2) ID_best_model.R - Identify the model with the largest log marginal likelihood from the stepping stone sampler files. Saves model choice as modelchoice.txt in the `home/data/BAGERresults/$tissuecode` directory

```
R --vanilla < ID_best_model.R
```


#### 3) Create.BAGER.Summary.File.R - Collect BAGER summary statistics (foldSD change between ancestral and descendent reconstructions, medians of both reconstructions, and the Bayesian Posterior Probability of Divergence/BPPD [see script for description]) for every branch into one Summary file, called `home/data/BAGERresults/$tissuecode/$DATEBAGERSummary.txt`.

```
R --vanilla < Create.BAGER.Summary.File.R 
```

#### 4) IdentifyShiftsInAllLineages.Rmd - Analyzes Ancestral Transcriptome Reconstructions to identify genes with expression shifts across the phylogeny. Also, finds enriched pathways in each lineage and visualizes transcriptome evolution.


Best used in R studio.


## Simulate gene expression evolution across your tree. This portion is not yet finished.

**Using the `SimulateAGER.Rmd` file, you may simulate data that matches your data's expression patterns.**


## Troubleshooting on test example files

To run these scripts on the example test files, copy the files in scripts/test to the data directory and follow the above steps.

### written by Erin Fry
### Last modified: August 21 2017