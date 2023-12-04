# SPA<sub>GRM</sub> introduction

SPA<sub>GRM</sub> is a scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies (GWAS). In the paper ```A scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank (to be updated)```, we applied SPA<sub>GRM</sub> to analyze 79 longitudinal traits extracted from UK Biobank primary care data. As a universal analysis framework, we also evaluated SPA<sub>GRM</sub>'s performance in quantitative and binary trait analysis. 

SPA<sub>GRM</sub> is a two-step method to control for sample relatedness in large-scale cohort like many other popular methods, such as [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html), [SAIGE](https://saigegit.github.io/SAIGE-doc/), [fastGWA](https://yanglab.westlake.edu.cn/software/gcta/#Overview), [REGENIE](https://rgcgithub.github.io/regenie/), [GATE](https://github.com/weizhou0/GATE), and [POLMM](https://github.com/WenjianBI/POLMM). In step 1, SPA<sub>GRM</sub> fits a null model to adjust for the effect of covariates on phenotypes and calculate model residuals. It is optional, rather than required, to incorporate a random effect into null model fitting to characterize the sample relatedness. In step 2, SPA<sub>GRM</sub> associates the trait of interest to a single genetic variant and obtain GWAS results by a retrospective strategy: that is treating the genotypes as random variables. Saddlepoint approximation is applied in SPA<sub>GRM</sub> that can greatly increase the accuracy to analyze low-frequency and rare variants, especially if the phenotypic distribution is unbalanced. To avoid redundant computations in step 2, SPAGRM uses the genetic relationship matrix (GRM) and identity by descent (IBD) probabilities to approximate the discrete joint distribution of genotype vectors in advance (illustrated as step 0). See ```A scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank (to be updated)``` for more details about workflow of SPA<sub>GRM</sub>.

![plot](https://github.com/Fantasy-XuHe/SPAGRM/blob/main/pictures/workfolw%20of%20SPAGRM.png)

SPA<sub>GRM</sub> is now implemented in the [GRAB package](https://wenjianbi.github.io/grab.github.io/)

Detailed documentation is available at [SPAGRM online tutorial](https://fantasy-xuhe.github.io/SPAGRM.github.io/)

# SPA<sub>GRM</sub> reproducibility

Scripts to reproduce the experiments performed for the SPA<sub>GRM</sub> manuscript:

A scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank.

## Simulation studies

## Real data analysis
