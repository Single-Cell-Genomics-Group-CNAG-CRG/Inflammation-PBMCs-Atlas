#!/bin/bash

parallel --delay 2 -j 1 eval "sbatch --qos=normal --time=11:00:00 --job-name='{1}' --mem=500G -c 12 --output {1}.log --wrap=\'papermill /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/01_run_ULMsignatures_SPECTRAfactors.ipynb /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/executed_notebooks/01_run_ULMsignatures_SPECTRAfactors__{1}.ipynb -k python3 -p celltype {1}\'" ::: T_CD8_NonNaive T_CD8_Naive Mono T_CD4_Naive T_CD4_NonNaive ILC B DC pDC UTC Plasma


## Subsetting donors
parallel --delay 2 -j 1 eval "sbatch --qos=normal --time=06:00:00 --job-name='{1}' --mem=500G -c 12 --output {1}.log --wrap=\'papermill /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/01_run_MLMsignatures_SPECTRAfactors_noWeights_subsetPatients.ipynb /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/executed_notebooks/01_run_MLMsignatures_SPECTRAfactors_noWeights_subsetPatients__{1}.ipynb -k python3 -p celltype {1}\'" ::: T_CD8_NonNaive T_CD8_Naive Mono T_CD4_Naive T_CD4_NonNaive ILC B DC pDC UTC Plasma 


# NOT CHECKING Progenitors Cycling_cells

# Pseudobulk before UML :)
# Corrected

parallel --delay 2 -j 1 eval "sbatch --qos=normal --time=11:00:00 --job-name='{1}' --mem=400G -c 12 --output {1}.log --wrap=\'papermill /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/01_run_ULMsignaturesL1_Corr_SPECTRAfactors.ipynb /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/executed_notebooks/01_run_ULMsignaturesL1_Corr_SPECTRAfactors__{1}.ipynb -k python3 -p celltype {1}\'" ::: T_CD8_NonNaive T_CD8_Naive Mono T_CD4_Naive T_CD4_NonNaive ILC B DC pDC UTC Plasma

parallel --delay 2 -j 1 eval "sbatch --qos=normal --time=11:00:00 --job-name='{1}' --mem=400G -c 12 --output {1}.log --wrap=\'papermill /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/01_run_ULMsignaturesL2_Corr_SPECTRAfactors.ipynb /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/executed_notebooks/01_run_ULMsignaturesL2_Corr_SPECTRAfactors__{1}.ipynb -k python3 -p celltype {1}\'" ::: T_CD8_NonNaive Mono T_CD4_NonNaive ILC B DC UTC Plasma


# UnCorrected
parallel --delay 2 -j 1 eval "sbatch --qos=normal --time=11:00:00 --job-name='{1}' --mem=400G -c 12 --output {1}.log --wrap=\'papermill /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/01_run_ULMsignaturesL1_UnCorr_SPECTRAfactors.ipynb /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/executed_notebooks/01_run_ULMsignaturesL1_UnCorr_SPECTRAfactors__{1}.ipynb -k python3 -p celltype {1}\'" ::: T_CD8_NonNaive T_CD8_Naive Mono T_CD4_Naive T_CD4_NonNaive ILC B DC pDC UTC Plasma

parallel --delay 2 -j 1 eval "sbatch --qos=normal --time=11:00:00 --job-name='{1}' --mem=400G -c 12 --output {1}.log --wrap=\'papermill /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/01_run_ULMsignaturesL2_UnCorr_SPECTRAfactors.ipynb /scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/06_inflammation_signatures/executed_notebooks/01_run_ULMsignaturesL2_UnCorr_SPECTRAfactors__{1}.ipynb -k python3 -p celltype {1}\'" ::: T_CD8_NonNaive Mono T_CD4_NonNaive ILC B DC UTC Plasma