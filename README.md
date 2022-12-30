# IR - BW project (add title)

## GWAS summary statistics used:

In this study we explored the genetic association between insulin secretion and insulin resistance (exposures) and birthweight (outcome). We performed four analysis that can be divided into two groups depending on the type of data that we utilized:

a) analysis that use selected instrumented variants: genetic risk scores (GRS) and common mendelian randomization analysis with 2SMR package.
b) analysis that use all genome-wide association studies: genetic correlations with High Definition Likelihood (HDL) and Linkage Disequilibrium Score (LDSC) and CAUSE MR package.

For each exposure we will here report which data we used in each analysis:

### Exposure: insulin secretion:

For insulin sensitivity we used three traits: 

a) Corrected insulin response (CIR)
b) Corrected insulin response adjusted for insulin sensitivity index (CIRadjISI)
c) Deposition Index (DI)

All GWAS summary statitics are from Prokopenko et al can be found in MAGIC consortium (https://magicinvestigators.org/downloads/).

Importantly, for GRS and common MR analysis with 2SMR packages we used the variants associated with insulin secretion reported in Table 1 from Prokopenko's paper:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3974640/table/pgen-1004235-t001/?report=objectonly. This data was used for every traits with the logic that these 8 variants are the instrumental variants for the strongest insulin secretion traits (CIR), so we also used them for CIRadjISI and DI analysis to see if the results found with CIR are reproducible with similar traits. Go to Supplementary Note in our publication for a lengthy explanation and the sensitivity analysis we performed regarding this issue.

### Exposure: insulin resistance

For insulin resistance we also used three traits:

a) fasting insulin (FI)
b) fasting insulin adjusted for body mass index (FIadjBMI)
c) Insulin Resistance (IR)

For analysis that used all GWAS summary statistics (genetic correlations and CAUSE MR) we used the whole summary statistics for FI, FIadjBMI and HOMA-IR as a proxy for IR. FI data can be downladed here: https://magicinvestigators.org/downloads/FI_combined_1000G_density_formatted_21-03-29.txt.gz; FIadjBMI can be downloaded here: https://magicinvestigators.org/downloads/MAGIC1000G_FI_EUR.tsv.gz; and HOMA-IR data can be downloaded here: https://magicinvestigators.org/downloads/MAGIC_ln_HOMA-IR.txt.

For analysis that used instrumented variants we used the variants reported in the tables of FI and FIadjBMI publications -Lagou et al and Che et al-. For IR we decided to utilized the 53 IR variants identified in Lotta et al 2016 using an integrative genomic approach. (See Supplementary Note in publication and curating_data folder for more information).

### Outcome: birthweight 

For birthweight we used 5 traits coming from 2 different sources. Since this data was used as outcome, we just queried the summary statistics needed for each exposure variant and the birthweight GWAS.

We used 3 unadjusted traits from Juliusdottir et al:

a) fetal effect on birthweight (EGG+UKBB+Icelandic GWAS meta-analysis)
b) maternal effects on birthweight (EGG+UKBB+Icelandic GWAS meta-analysis)
c) paternal effects on birthweight (only Icelandic GWAS available)

Additionally, we used 2 adjusted traits from Warrington et al to disentangle the associations that are driven by the fetus' or by the mother's genetics

a) fetal effects adjusted for maternal genotype
b) maternal effects adjsted for fetal genotype

Juliusdottir data can be obtained here: https://www.decode.com/summarydata/

Warrington data can be obtained here: http://egg-consortium.org/birth-weight-2019.html

## Folder structure:


