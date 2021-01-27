
# MTXmodel User Manual #

This R package is built for metatranscriptomics (MTX) modeling based on [MaAsLin2](http://huttenhower.sph.harvard.edu/maaslin2). It integrates feature-specific covariates to determine multivariable association between metadata and microbial MTX features. MTX abundance changes are highly affected by underlying differences in metagenomic abundances (i.e. gene copy number). This package can adjust the DNA abundance as a continuous covariate for a given feature in the models for differential expression analysis in microbial communities. 


If you use the MTXmodel package, please cite our manuscript:
Zhang et al. (2021+). "Statistical approaches for differential expression analysis in metatranscriptomics" (In Submission).

And feel free to link it to your Methods: http://huttenhower.sph.harvard.edu/mtx2021

--------------------------------------------

## Contents ##
* [Description](#description)
* [Requirements](#requirements)
* [Installation](#installation)
* [How to Run](#how-to-run)
    * [Input Files](#input-files)
    * [Output Files](#output-files)
    * [Run a Demo](#run-a-demo)
    * [Options](#options)

## Description ##

MTXmodel keeps all the features and functions of [MaAsLin2](https://github.com/biobakery/Maaslin2) and additionally add new functions which can be used for adjusting feature-specific covariates.

## Requirements ##

MTXmodel is an R package that can be run on the command line or as an R function.

## Installation ##

MTXmodel can be run from the command line or as an R function.

If only running from the command line, you do not need to install the MTXmodel package but you will need to install the MTXmodel dependencies.

### From command line ###

1. Download the source: [MTXmodel.master.zip](https://github.com/biobakery/mtx2021/archive/master.zip)
2. Decompress the download: 
    * ``$ tar xzvf MTXmodel-master.zip``
3. Install the Bioconductor dependencies edgeR and metagenomeSeq. 
4. Install the CRAN dependencies:
    * ``$ R -q -e "install.packages(c('lmerTest','pbapply','car','dplyr','vegan','chemometrics','ggplot2','pheatmap','hash','logging','data.table','MuMIn','glmmTMB','MASS','cplm','pscl'), repos='http://cran.r-project.org')"``
5. Install the MTXmodel package (only required if running as an R function): 
    * ``$ R CMD INSTALL MTXmodel-master``


## How to Run ##

MTXmodel can be run from the command line or as an R function. Both methods require the same arguments, have the same options, and use the same default settings.

### Input Files ###

MTXmodel requires three input files.

1. Data (or features) file
    * This file is tab-delimited.
    * Formatted with features as columns and samples as rows.
    * The transpose of this format is also okay.
    * Possible features in this file include taxonomy or genes or pathways.
2. Metadata file
    * This file is tab-delimited.
    * Formatted with features as columns and samples as rows.
    * The transpose of this format is also okay.
    * Possible metadata in this file include gender or age.
3. Covariate data of features file
    * This file is tab-delimited.
    * Formatted with features as columns and samples as rows.
    * The transpose of this format is also okay.
    * Possible data in this file include DNA abundance of genes or pathways.

The data file can contain samples not included in the metadata file
(along with the reverse case). For both cases, those samples are not 
included in both files will be removed from the analysis. 
Also the samples do not need to be in the same order in the two files.

NOTE: If running MTXmodel as a function, the data and metadata 
inputs can be of type ``data.frame`` instead of a path to a file.

### Output Files ###

MTXmodel generates the same types of output files with MaAsLin2: data and visualization. See more details in [MaAsLin2 manual](https://github.com/biobakery/Maaslin2#output-files) 

### Run a Demo ###

Example input files can be found in the ``inst/extdata`` folder 
of the MTXmodel source. The files provided were generated from
the HMP2 data which can be downloaded from https://ibdmdb.org/ .

``HMP2_pwyRNA.tsv``: is a tab-delimited file with pathways as columns and samples as rows. It is a subset of the pathway file so it just includes the pathway RNA abundances for all samples.

``HMP2_pwyDNA.tsv``: is a tab-delimited file with pathways as columns and samples as rows. It is a subset of the pathway file so it just includes the pathway DNA abundances for all samples.

``HMP2_metadata.tsv``: is a tab-delimited file with samples as rows and metadata as columns. It is a subset of the metadata file so that it just includes some of the fields.


#### Command line ####

``$ MTXmodel.R ../inst/extdata/HMP2_pwyRNA.tsv ../inst/extdata/HMP2_metadata.tsv demo_output --min_abundance 0 --min_prevalence 0.0 --max_significance 0.25 --min_variance 0.0 --correction BH --standardize TRUE --normalization NONE --transform LOG --analysis_method LM --cores 1 --fixed_effects diagnosis,dysbiosisCD,dysbiosisUC,dysbiosisnonIBD,antibiotics,age --random_effects subject,site --plot_heatmap FALSE --plot_scatter FALSE --reference diagnosis,nonIBD --rna_dna_flt local --input_dnadata  ../inst/extdata/HMP2_pwyDNA.tsv
``

* Make sure to provide the full path to the MTXmodel executable (i.e., ./R/MTXmodel.R).
* In the demo command:
    * ``HMP2_pwyRNA.tsv`` is the path to your data (or features) file
    * ``HMP2_pwyDNA.tsv`` is the path to your paired DNA data of features file
    * ``HMP2_metadata.tsv`` is the path to your metadata file
    * ``demo_output`` is the path to the folder to write the output


#### In R ####

```{r}
library(MTXmodel)
input_data <- system.file(
    'extdata','HMP2_pwyDNA.tsv', package="MTXmodel")
input_metadata <-system.file(
    'extdata','HMP2_metadata.tsv', package="MTXmodel")
input_dnadata <- system.file(
    'extdata','HMP2_pwyRNA.tsv', package="MTXmodel")
fit_data <- MTXmodel(
    input_data, input_metadata, 'demo_output', transform = "LOG",
    fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
    random_effects = c('site', 'subject'),
    reference = "diagnosis,nonIBD",
    normalization = 'NONE',
    standardize = FALSE,
    input_dnadata = input_dnadata
    )
```

##### Session Info #####

Session info from running the demo in R can be displayed with the following command.

```{r}
sessionInfo()
```

### Options ###

Run MTXmodel help to print a list of the options and the default settings.


$ MTXmodel.R --help

Usage: ./R/MTXmodel.R [options] <data.tsv> <metadata.tsv> <output_folder>


Options:

	-h, --help
		Show this help message and exit
		
	-a MIN_ABUNDANCE, --min_abundance=MIN_ABUNDANCE
		The minimum abundance for each feature [ Default: 0 ]

	-p MIN_PREVALENCE, --min_prevalence=MIN_PREVALENCE
		The minimum percent of samples for which a feature is detected at minimum abundance [ Default: 0.1 ]

	-b MIN_VARIANCE, --min_variance=MIN_VARIANCE
		Keep features with variances greater than value [ Default: 0 ]

	-s MAX_SIGNIFICANCE, --max_significance=MAX_SIGNIFICANCE
		The q-value threshold for significance [ Default: 0.25 ]

	-n NORMALIZATION, --normalization=NORMALIZATION
		The normalization method to apply  [ Default: TSS ] [ Choices: TSS, CLR, CSS, NONE, TMM ]

	-t TRANSFORM, --transform=TRANSFORM
		The transform to apply [ Default: LOG ] [ Choices: LOG, LOGIT, AST, NONE, PA ]

	-m ANALYSIS_METHOD, --analysis_method=ANALYSIS_METHOD
		The analysis method to apply [ Default: LM ] [ Choices: LM, CPLM, NEGBIN, ZINB, LOGIT ]

	-r RANDOM_EFFECTS, --random_effects=RANDOM_EFFECTS
		The random effects for the model,  comma-delimited for multiple effects  [ Default: none ]

	-f FIXED_EFFECTS, --fixed_effects=FIXED_EFFECTS
		The fixed effects for the model,  comma-delimited for multiple effects  [ Default: all ]

	-c CORRECTION, --correction=CORRECTION
		The correction method for computing  the q-value [ Default: BH ]

	-z STANDARDIZE, --standardize=STANDARDIZE
		Apply z-score so continuous metadata are on  the same scale [ Default: TRUE ]

	-l PLOT_HEATMAP, --plot_heatmap=PLOT_HEATMAP
		Generate a heatmap for the significant  associations [ Default: TRUE ]

	-i HEATMAP_FIRST_N, --heatmap_first_n=HEATMAP_FIRST_N
		In heatmap, plot top N features with significant  associations [ Default: 50 ]

	-o PLOT_SCATTER, --plot_scatter=PLOT_SCATTER
		Generate scatter plots for the significant  associations [ Default: TRUE ]

	-e CORES, --cores=CORES
		The number of R processes to  run in parallel [ Default: 1 ]

	-d REFERENCE, --reference=REFERENCE
		The factor to use as a reference for a variable with more than two levels provided as a string of 'variable,reference' semi-colon delimited for multiple variables [ Default: NA ]

	-x INPUT_DNADATA, --input_dnadata=INPUT_DNADATA
		The DNA abundance for each feature [ Default: none ]

	-y RNA_DNA_FLT, --rna_dna_flt=RNA_DNA_FLT
		Filtering features/samples based on the detectable abundance of RNA and covariate DNA per feature [ Default: global ]
		[ Choices: none, global, local, strict]
		none: do not apply filtering
		global: ignore features that are not detected at both DNA and RNA levels across all samples
		local: ignore the sample where both feature's DAN and feature's RNA are not detected
		strict: ignore the sample where either feature's DAN or feature's RNA is not detected
