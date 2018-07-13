# PINES
Phenotype-informed tissue weighting improves prediction of disease-relevant noncoding variants

Functional characterization of the noncoding genome is essential for the biological understanding of gene regulation and disease. Here, we introduce the computational framework PINES (Phenotype-Informed Noncoding Element Scoring) which predicts the functional impact of noncoding variants by integrating epigenetic annotations in a phenotype-dependent manner. A unique feature of PINES is that analyses may be customized towards genomic annotations from cell types of the highest relevance given the phenotype of interest.

Requirements: R (including library optparse), tabix

The epigenetic annotations needed to run PINES can be found at: https://zenodo.org/record/1311246#.W0h6y9hKiLI. 

The script PINES.R accepts the following command-line arguments:

--f: path to the text file containing the RS IDs (one per line) of the variants to be scored
--e: path to the text file containing the RS IDs (one per line) of the variants used for enrichment-based weighting of the epigenetic annotations
--w: manual weighting constant
--t: tissue to upweight: select from neuro, liver, immune, or heart
--c: number of cores to use for parallel computations
--p: path to the annotation folder

Examples:

1. PINES run by specifying a set of variants for enrichment-based weighting of the epigenetic annotations: 
Rscript PINES.R --p /path/to/annotations_folder/ --f /path/to/test_variants.txt --e /path/to/GWAS_variants.txt  --c 2

2. PINES run by specifying a pre-defined set of cell types and a manual weighting constant:
Rscript PINES.R --p /path/to/annotations_folder/ --f /path/to/test_variants.txt --t immune  --w 4 --c 2

The script draw_enrichment_heatmap.R plots the annotation weight matrix that is generated by PINES.R if a set of variants was specified for enrichment-based weighting. The script draw_enrichment_heatmap.R accepts the following command-line arguments:

--f: path to the PINES output file copy_of_feature_enrichment_pvalues.RData
--p: path to the annotation folder

Example:

Rscript draw_enrichment_heatmap.R --p /path/to/annotations_folder/ --f /path/to/copy_of_feature_enrichment_pvalues.RData

PINES can also be run through an interactive web portal at: http://genetics.bwh.harvard.edu/pines/
