# PINES
Phenotype-informed tissue weighting improves prediction of disease-relevant noncoding variants

Functional characterization of the noncoding genome is essential for the biological understanding of gene regulation and disease. Here, we introduce the computational framework PINES (Phenotype-Informed Noncoding Element Scoring) which predicts the functional impact of noncoding variants by integrating epigenetic annotations in a phenotype-dependent manner. A unique feature of PINES is that analyses may be customized towards genomic annotations from cell types of the highest relevance given the phenotype of interest.

Requirements: R, tabix

The epigenetic annotations needed to run PINES can be found at: https://zenodo.org/record/1228512#.WuBft1MvyLI. The path to the downloaded and unzipped annotation folder needs to be specified at the beginning of the R script under the "annotation_path" variable.

PINES can also be run through an interactive web portal at: http://genetics.bwh.harvard.edu/pines/
