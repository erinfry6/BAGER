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
					
						$tissuecode
				
				home/code
				
				home/BayesTraits
				
				
Place the contents of this repository's `code` directory and `analysis/4_IdentifyShiftsInAllLineages.Rmd` in the `code` directory. Place data to be analyzed in `data/forBAGER` and create directories for each tissue in `data/BAGERresults`.

It may be useful to consult the [BayesTraits Manual](http://www.evolution.rdg.ac.uk/BayesTraitsV3/Files/BayesTraitsV3.Manual.pdf) to understand the phylogenetic models used to reconstruct ancestral traits.


## Input files formats

The following input files are required for this analysis and should be placed in `data/forBAGER`. (Examples of these formats can be found in `demo` in part created by `code/creating_demo.R`):

 - tissuecode_exp.txt: a tab delimited file containing expression data formatted according to the BayesTraits Manual: Expression data file with no column names. First column: names of the samples that coordinate with the Nexus tree file. All subsequent columns are the expression data for each gene. *Be sure to keep another document (_genesincluded) that documents the gene in each column and which gene number corresponds to which gene.*

 - tissuecode_tree.tree: an ultrametric nexus formatted phylogeny formatted according to the BayesTraits Manual.
 
 - tissuecode_genesincluded.txt: a tab delimited file containing expression data identical to that in tissuecode_exp.txt, but formatted differently. This file should have column and row names. Columns are the samples and rows are each gene, labeled with its EnsemblID.


**tissuecode is the name of the tissue you are analyzing that will be an input for the scripts below**


## Modify the scripts

 - **Please see the instructions at the top of each individual script for modification requirements.**

 - If you would like to modify the models tested by BayesTraits from the default models in this pipeline, you will need to specify which models in `create_model_files.sh` , `identify_best_model.sh` , `ID_best_model.R`, and `Create.BAGER.Summary.File.R`.

 
## Run BAGER


#### 1) 1_ancestral_reconstruction.sh - Reconstruct transcriptomes under each model. The default recommended models are coded in, but if you would like to change or add models, you will have to do so.
I typically run the code one tissue at a time, broken into 4 or 5 chunks specified by begingene and endgene numbers to decrease run time per tissue. Begingene and endgene are integers of the beginning and start gene number you would like to reconstruct.

```
./1_ancestral_reconstruction.sh tissuecode begingene endgene
```


#### 2) 2_ID_best_model.R - Identify the model with the largest log marginal likelihood from the stepping stone sampler files. Saves model choice as modelchoice.txt in the `home/data/BAGERresults/$tissuecode` directory

```
R --vanilla < 2_ID_best_model.R
```


#### 3) 3_Create.BAGER.Summary.File.R - Collect BAGER summary statistics (foldSD change between ancestral and descendent reconstructions, medians of both reconstructions, and the Bayesian Posterior Probability of Divergence/BPPD [see script for description]) for every branch into one Summary file, called `home/data/BAGERresults/$tissuecode/$DATEBAGERSummary.txt`.

```
R --vanilla < 3_Create.BAGER.Summary.File.R 
```

#### 4) 4_IdentifyShiftsInAllLineages.Rmd - Analyzes Ancestral Transcriptome Reconstructions to identify genes with expression shifts across the phylogeny. Also, finds enriched pathways in each lineage and visualizes transcriptome evolution. You will need to run step 5 before knitting this file.


Best used in R studio.

#### 5) 5_consolidate_enrichments.sh - Concatenates all significant enrichments from enrichment outputs in previous step into one file in `home/data/BAGERresults/$tissuecode/Enrichments/enrichmentresults.txt`

```
./5_consolidate_enrichments.sh $tissuecode
```


## Troubleshooting on demo example files

To run these scripts on the example demo files provided, copy the files in `demo` to `data/forBAGER` and follow the above steps.


## Simulate gene expression evolution across your tree. This portion is not yet finished.

**Using the `SimulateAGER.Rmd` file, you may simulate data that matches your data's expression patterns.**


### written by Erin Fry
### Last modified: August 23 2017