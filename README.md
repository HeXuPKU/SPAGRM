# SPA<sub>GRM</sub> introduction

SPA<sub>GRM</sub> is a scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies (GWAS). In the paper ```A scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank (to be updated)```, we applied SPA<sub>GRM</sub> to analyze 79 longitudinal traits extracted from UK Biobank primary care data. As a universal analysis framework, we also evaluated SPA<sub>GRM</sub>'s performance in quantitative and binary trait analysis. 

**SPA<sub>GRM</sub> is now implemented in the [GRAB package](https://wenjianbi.github.io/grab.github.io/). Please click here to download.**

**Detailed documentation about how to use SPA<sub>GRM</sub> is available at [SPA<sub>GRM</sub> documentation](https://fantasy-xuhe.github.io/SPAGRM.github.io/).**

**Summary statistics of 79 longitudinal traits extracted from UK Biobank primary care data is available at [here](https://zenodo.org/records/10242062).**

# About this repository

This repository contains: 
1) old version of SPA<sub>GRM</sub>;
2) figures of GWAS results of 79 longitudinal traits;
3) materials for reproducing the experiments.

# SPA<sub>GRM</sub> reproducibility

Scripts to reproduce the experiments performed for the SPA<sub>GRM</sub> manuscript:

```A scalable, accurate, and universal analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank (to be updated)```

## Real data analysis
1. extract longitudinal traits from UK Biobank primary care data
```
In real_data\1.extract_pheno:
# phecode-2023-05-10XH.R defines the Read v2 and CTV3 code terms for 79 longitudinal traits.
# real_data_extra-2023-05-01XH.R preprocesses the gp_clinical table.
```
2. preprocess longitudinal traits
```
```
3. fit the null model via WiSER package
```
```
4. GWAS analysis
- SPA<sub>GRM</sub> analysis 
```
```
- TrajGWAS analysis
```
```
5. PRS adjustment (optional)
```
```

## Simulation study
