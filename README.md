# **Blood DNA Methylation Signature for Incident Dementia: Evidence from Longitudinal Cohorts** 
Wei Zhang, Juan I. Young, Lissette Gomez, Michael A. Schmidt, David Lukacsovich, Achintya Varma, Brian Kunkle, X. Steven Chen, Eden R. Martin, Lily Wang

## Discription

Alzheimer’s disease and related dementias pose a significant public health challenge, especially as the population ages. Dementia cases are often underreported, highlighting the need to identify individuals at risk early. However, distinguishing between molecular changes that precede dementia onset and those resulting from the disease is challenging with cross-sectional studies. To address this, we studied blood DNA methylation (DNAm) differences and incident dementia in two large longitudinal cohorts: the Offspring cohort of the Framingham Heart Study (FHS) and the Alzheimer’s Disease Neuroimaging Initiative (ADNI) study.

We analyzed blood DNAm samples from over 1,000 cognitively unimpaired subjects. FHS participants (n = 907) were followed for up to 7.72 years after blood sample collection at Exam 9; ADNI participants (n = 216) were followed for up to 11.11 years after their initial visits. Mean ages at sample collection were 72.03 years in FHS and 76.73 years in ADNI. Meta-analysis of results from cox regression models identified 44 CpGs and 44 differentially methylated regions consistently associated with time to dementia in both cohorts. Our integrative analysis identified early processes in dementia, such as immune responses and metabolic dysfunction. External validations with two independent datasets, the Australian Imaging, Biomarkers, and Lifestyle (AIBL) study and the AddNeuroMed study, showed significant discriminatory classification of dementia samples versus controls using methylation risk scores based on 44 dementia-associated CpGs. These findings demonstrate that DNA methylation offers a promising pathway for early detection and prevention of dementia in at-risk populations. 

### 1. Preprocessing of DNA methylation data

DNA methylation samples from both FHS and ADNI were measured using the same Illumina HumanMethylation EPIC v1 bead chips. Quality control of both probes and samples were performed. 

| File and folder                                              | Description                                       |
| ------------------------------------------------------------ | ------------------------------------------------- |
| [code/ADNI/preprocessing](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/ADNI/preprocessing) | Code for preprocessing of ADNI DNAm blood samples |
| [code/Framingham/preprocessing](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/Framingham/preprocessing) | Code for preprocessing of FHS9 DNAm blood samples |
### 2. Association of DNA methylation at individual CpGs with dementia

To evaluate the relationship between incident dementia and DNA methylation, we conducted Cox proportional regression analyses on both FHS and ADNI datasets separately, via the coxph function in the survival R package.  

| File and folder                                              | Description                           |
| ------------------------------------------------------------ | ------------------------------------- |
| [code/ADNI/cpg_test/cox_cpg_test_ADNI.Rmd](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/ADNI/cpg_test/cox_cpg_test_ADNI.Rmd) | Code for individual CpGs test of ADNI |
| [code/Framingham/cpg_test/cox_cpg_test_FHS9.Rmd](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/Framingham/cpg_test/cox_cpg_test_FHS9.Rmd) | Code for individual CpGs test of FHS9 |

### 3. Meta-analysis

To meta-analyze individual CpG results across both the FHS and ADNI datasets, we used the inverse-variance weighted fixed-effects model, implemented in the meta R package.

| File and folder                                              | Description            |
| ------------------------------------------------------------ | ---------------------- |
| [code/meta_analysis/meta_analysis.Rmd](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/meta_analysis/meta_analysis.Rmd) | Code for meta-analysis |
### 4. Pathway analysis

To identify biological pathways enriched with significant DNA methylation differences, we used the methylRRA function in the methylGSA R package[^1].

| File and folder                                              | Description               |
| ------------------------------------------------------------ | ------------------------- |
| [code/pathway_analysis](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/pathway_analysis) | Code for pathway analysis |

### 5. Integrative analyses with gene expression, genetic variants, and brain-to-blood correlations

To evaluate the effect of DNA methylation on the expression of nearby genes, we overlapped our dementia-associated CpGs, including both significant individual CpGs and those located within DMRs, with eQTm analysis results in Supplementary Tables 2 and 3 of Yao et al (2021)[^2].

To assess the correlation of dementia-associated CpGs and DMRs methylation levels in blood and brain samples, we used the London cohort, which consisted of 69 samples with matched PFC and blood samples.

| File and folder                                              | Description                          |
| ------------------------------------------------------------ | ------------------------------------ |
| [code/check_overlap](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/check_overlap) | Code for integrative analyses        |
| [code/brain_blood_corr](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/brain_blood_corr) | Code for brain-to-blood correlations |

### 6. Sensitivity analysis

**Covariate analysis**: We evaluated if dementia risk factors would likely confound the DNA methylation to dementia associations we observed. To this end, we first performed regression analysis to assess the association between dementia-associated CpGs and dementia risk factors collected by the Framingham study.

**Family structure analysis**: To evaluate the impact of family structure in the discovery of significant CpGs, ww compared the results a model that accounts for family relationships in the Cox regression model using a kinship matrix. 

| File and folder                                              | Description                        |
| ------------------------------------------------------------ | ---------------------------------- |
| [code/covariates_analysis](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/covariates_analysis) | Code for covariate analysis        |
| [code/Framingham/cpg_test/kinship.Rmd](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/Framingham/cpg_test/kinship.Rmd) | Code for family structure analysis |

### 7. External validation

The external validation analysis of the methylation risk scores (MRS) included 491 whole blood DNA methylation samples from the AIBL study and 171 whole blood DNA methylation samples from the cross-European AddNeuromed study, which were downloaded from the GEO database (accession: GSE153712 and GSE144858). The preprocessing of DNA methylation samples from AIBL and AddNeuromed were previously described in Silva et al. (2022)[^3].

| File and folder                                              | Description           |
| ------------------------------------------------------------ | --------------------- |
| [code/MRS_analysis](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/MRS_analysis) | Code for MRS analysis |

## For reproducible research

To perform the analysis, begin by installing the packages found in `session_info.R` ([Link to the script](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/session_info.R)). Then, load the auxiliary functions from folder `code/utility` ([Link to the file](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/utility)). Follow the sequence provided in the Description to conduct the analysis.

## Acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

## Reference

[^1]: Ren, X. & Kuan, P.F. methylGSA: a Bioconductor package and Shiny app for DNA methylation data length bias adjustment in gene set testing. *Bioinformatics* **35**, 1958-1959 (2019)
[^2]: Yao, C. *et al.* Epigenome-wide association study of whole blood gene expression in Framingham Heart Study participants provides molecular insight into the potential role of CHRNA5 in cigarette smoking-related lung diseases. *Clin Epigenetics* **13**, 60 (2021)  
[^3]: T, C.S. *et al.* Cross-tissue analysis of blood and brain epigenome-wide association studies in Alzheimer's disease. *Nat Commun* **13**, 4852 (2022)

