# Pan-Cancer Spatial Atlas of TLS

This repository includes the codes used for generating the main results presented in the paper Pan-Cancer Spatial Atlas of Tertiary Lymphoid Structures. The repository includes data preprocessing, statistical analysis, visualization, and supporting utilities needed to reproduce the study outputs. Please refer to the folder code/ for related code. Here we provide the outline of analysis covered by each script:

`Fig1.R` covers the following analyses:  
> - Pan-cancer quantification of TLS-positive versus TLS-negative samples.  
> - Comparison of TLS density and TLS count distributions across cancer types.  
> - Characterization of TLS spatial locations and enrichment patterns by cancer type.

`Fig2.R` covers the following analyses:  
> - Composition of TLS maturation states.  
> - Comparative analysis of TLS score and marker expression across maturation states.  
> - Association of maturation states with spatial location, TLS density, and cancer type.

`Fig3.R` covers the following analyses:  
> - Marker-level gene expression heatmaps across TLS maturation states.  
> - Comparison of transcriptional programs linked to TLS maturation.

`Fig4.R` covers the following analyses:  
> - Spatial gradient analysis using distance from TLS as a continuous spatial ordering metric.  
> - Identification of significantly dynamic pathways/genes along the TLS-distance gradient and trend-direction classification.  
> - Cross-cancer summarization and visualization of conserved TLS-distance-associated pathway dynamics.

`Fig5.R` covers the following analyses:  
> - Cohort-specific TLS maturation profiling.  
> - Integration of TLS composition with clinical and molecular features.  
> - Spatial context analysis of maturation states by tissue region.

`Fig6.R` covers the following analyses:  
> - TLS maturation composition analysis in H&E cohorts.  
> - Construction of TLS composition-based groups and evaluation of prognostic associations across cancers.  
> - Cohort-specific treatment/outcome analyses comparing TLS composition by response.

`Fig1.sh` covers the following analysis workflow:  
> - TLS detection from ST samples.  
> - Marker scoring and phenotype inference to generate TLS-related embedding outputs.  
> - TLS clustering, segmentation, and extraction of TLS-level summary outputs.

`Fig6.sh` covers the following analysis workflow:  
> - TLS detection from H&E slides.  
> - TLS maturation classification.

For any questions, please leave your comment in GitHub or contact Kevin Cho (kcho2@mdanderson.org). We will help address the issues as soon as possible.

