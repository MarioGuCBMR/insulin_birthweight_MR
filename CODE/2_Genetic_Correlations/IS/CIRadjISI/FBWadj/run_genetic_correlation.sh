#! /bin/bash
#$ -S /bin/bash
#$ -N "CIRadjISI_FBWadj_GC"
#$ -cwd
#$ -pe smp 3
#$ -l mem_free=4G
#$ -l h_vmem=7G
#$ -l s_vmem=4G

conda activate r_env

Rscript /emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/2_Genetic_Correlations/0_installing_HDL_and_reference_data/HDL/HDL.run.R \
gwas1.df=/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/2_Genetic_Correlations/curated_files_4_gc/CIRadjISI_4_HDL.txt \
gwas2.df=/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/2_Genetic_Correlations/curated_files_4_gc/FBWadj_4_HDL.txt \
LD.path=/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/2_Genetic_Correlations/0_installing_HDL_and_reference_data/ref_panel/UKB_imputed_hapmap2_SVD_eigen99_extraction \
output.file=/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/2_Genetic_Correlations/output/CIRadjISI_FBWadj.Rout