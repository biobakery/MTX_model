MTXmodel.R ../inst/extdata/HMP2_pwyRNA.tsv ../inst/extdata/HMP2_metadata.tsv demo_output --min_abundance 0 --min_prevalence 0.0 --max_significance 0.25 --min_variance 0.0 --correction BH --standardize TRUE --normalization NONE --transform LOG --analysis_method LM --cores 1 --fixed_effects diagnosis,dysbiosisCD,dysbiosisUC,dysbiosisnonIBD,antibiotics,age --random_effects subject,site --plot_heatmap FALSE --plot_scatter FALSE --reference diagnosis,nonIBD --rna_dna_flt semi_strict --input_dnadata  ../inst/extdata/HMP2_pwyDNA.tsv

