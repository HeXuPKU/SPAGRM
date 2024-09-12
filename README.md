# SPA<sub>GRM</sub> introduction

SPA<sub>GRM</sub> is a scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies (GWAS). In the paper **A scalable and accurate analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in the UK Biobank (to be updated)**, we applied SPA<sub>GRM</sub> to analyze 79 longitudinal traits extracted from UK Biobank primary care data.

**SPA<sub>GRM</sub> is now implemented in the [GRAB package](https://wenjianbi.github.io/grab.github.io/). Please click here to download.**

**Detailed documentation about how to use SPA<sub>GRM</sub> is available at [SPA<sub>GRM</sub> documentation](https://hexupku.github.io/SPAGRM.github.io/).**

**Summary statistics of 79 longitudinal traits extracted from UK Biobank primary care data is available at [here](https://zenodo.org/records/10242062).**

# About this repository

This repository contains: 
1) old version of SPA<sub>GRM</sub>;
2) manhattan and QQ plots of GWAS results for 79 longitudinal traits;
3) materials for reproducing the experiments, including real data analyses and simulation studies.

# SPA<sub>GRM</sub> reproducibility

Scripts to reproduce the experiments performed for the SPA<sub>GRM</sub> manuscript:

**A scalable and accurate analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank (to be updated)**

## Real data analysis
1. extract longitudinal traits from UK Biobank primary care data
```
In real_data\1.extract_pheno:
# phecode-2023-05-10XH.R defines the Read v2 and CTV3 code terms for 79 longitudinal traits.
# real_data_extra-2023-05-01XH.R preprocesses the gp_clinical table.
```
2. preprocess longitudinal traits
```
In real_data\2.preprocess:
# XXX_preprocess.R extracts longitudinal traits from gp_clinical table and preprocesses each longitudinal traits.
```
3. fit the null model via WiSER package
```
In real_data\3.null_fitter:
# XXX_null_fitter.jl fits the null model and obtains model residuals for each longitudinal trait.
```
4. GWAS analysis
- SPA<sub>GRM</sub> analysis 
```
In real_data\4.1.SPAGRM:
# XXX_SPAGRM.R conducts the SPAGRM analysis with two parts: 1). pre-calculate of the joint distribution of genotypes and 2). conduct genome-wide association studies for each SNP.
```
- TrajGWAS analysis
```
In real_data\4.2.TrajGWAS:
# XXX_TrajGWAS.jl conducts the TrajGWAS analysis.
```
5. PRS adjustment (optional)
```
In real_data\5.PRS_adjustment (optional):
# 1.select_index_SNPs_from_results_of_SPAGRM.R selects independent SNPs that pass the given p value threshold based on the GWAS results of SPAGRM.
# 2.estimate_effect_sizes_XXX.jl estimate the effect sizes of index SNPs via WiSER package.
# 3.construct_LOCO-PRS.R constructs the Leave One Chromosome Out Polygenic Risk Scores (LOCO-PRS).
# 4.adjust_PRS_in_null_model_fitting_XXX.jl adds the PRS into the null model fitting.
# 5.second_round_of_GWAS_via_SPAGRM.R conducts the second round of GWAS via SPAGRM.
```

## Simulation study
```
In simulation:
# 1.simupheno contains R script to simulate longitudinal data.
# 2.typeIerror contains scripts to conduct 1e9 type I error simulations.
# 3.power contains scripts to conduct 1e3 power simulations.
```
