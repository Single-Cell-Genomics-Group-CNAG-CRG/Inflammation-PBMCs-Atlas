# *Interpretable Inflammation Landscape of Circulating Immune cells*
This is the Github repository including the code used in the paper entitled *Interpretable Inflammation Landscape of Circulating Immune cells*.

The manuscript is currently under a review process.

# Abstract
Inflammation is a biological phenomenon involved in a wide variety of physiological and pathological processes. Although a controlled inflammatory response is beneficial for restoring homeostasis, it can become unfavorable if dysregulated. In recent years, major progress has been made in characterizing acute and chronic inflammation in specific diseases. However, a global, holistic understanding of inflammation is still elusive. This is particularly intriguing, considering the crucial function of inflammation for human health and its potential for modern medicine if fully deciphered. Here, we leverage advances in the field of single-cell genomics to delineate the full spectrum of circulating immune cell activation underlying inflammatory processes during infection, immune-mediated inflammatory diseases and cancer. Our single-cell atlas of >6.5 million peripheral blood mononuclear cells from 1047 patients and 19 diseases allowed us to learn a comprehensive model of inflammation in circulating immune cells. The atlas expanded our current knowledge of the biology of inflammation of immune-mediated diseases, acute and chronic inflammatory diseases, infection and solid tumors, and laid the foundation to develop a precision medicine framework using unsupervised as well as explainable machine learning. Beyond a disease-centered analysis, we charted altered activity of inflammatory molecules in peripheral blood cells, depicting discriminative inflammation-related genes to further understand mechanisms of inflammation. Finally, we have laid the groundwork for developing precision medicine diagnostic tools for patients experiencing pathologic inflammation by learning a classifier for inflammatory diseases, presenting cells in circulation as a powerful resource for patient diagnosis.

![project_schema](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/Inflammation-PBMCs-Atlas/blob/main/graphical_abstract.png)

# Structure
The repository is organized into the following folder tree, which contains the precompiled notebooks that generated results shown in the manuscript.
```
├── 01_data_processing
|   └── SCGT00_CentralizedDataset
├── 02_cell_annotation
│   ├── 01_fromDatasets_to_CellLineages
│   ├── 02_fromCellLineages_to_CellTypes
│   │   ├── Step1
│   │   ├── Step2
│   │   ├── Step3
│   │   ├── Step4
│   │   └── Step5
│   ├── 03_characterizing_CellTypes
|   └── SCGT00_CentralizedDataset
│       ├── 01_fromDataset_to_CellLineages
│       └── 02_fromCellLineages_to_CellTypes
│           ├── Step1
│           ├── Step2
│           └── Step3
├── 03_downstream_analysis
│   ├── 01_compositional_analysis
│   ├── 02_gene_universe_definition
│   ├── 03_SEAcell_generation
│   ├── 04_integration_with_annotation
│   ├── 05_SPECTRA
│   ├── 06_inflammation_signatures
│   ├── 07_gene_regulatory_network
│   ├── 08_gene_importance
│   └── 09_patient_classifier
│       └── SCGT00_CentralizedDataset
├── 04_visualizing_final_embedding_space
│   └── SCGT00_CentralizedDataset
├── bin
└── external_reference_data

```

The repository is split into three main folders:

1. **Data processing**
Our single-cell atlas of >6.5 million peripheral blood mononuclear cells from 1047 patients and 19 diseases allowed us to learn a comprehensive model of inflammation in circulating immune cells. 
This folder contains all scripts necessary to perform the initial preprocessing of the dataset, which includes: Metadata homogeneization, Cell QC, Gene QC, dataset splitting into MAIN, VALIDATION and EXTERNAL and the initial integration of cells included in MAIN, using [scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html).
  * *SCGT00_CentralizedDataset*: From the unified data homogeneization, we focus on the "SCGT00" Centralized Dataset, and contains the same scripts as described before, to performn Cell / Gene QC and split the dataset into MAIN and EXTERNAL datasets, and the initial integration of cells included in the MAIN.

2. **Cell Annotation**:
This folder contains all scripts necessary to obtain the 64 different cell subpopulations described in the manuscript (Level 2). We employed a recursive top-down approach inspired by previous work done by La Manno et al. and Massoni-Badosa et al. Starting with more than 4M cells from MAIN dataset collected for the project, we divided the annotation into 2 major steps:
  * *01_fromDatasets_to_CellLineages*: We split our dataset into main Lineages.
    * Templates: To generate the scanpy object used as input for each annotation task, we executed the notebooks included in this directory.
  * *02_fromCellLineages_to_CellTypes*: We performed recursive steps to detected potential Doublets and Low quality cells, to be removed after validation. Apart from detecting all the celltypes of our interest, we also captured cells resembling Platelets or Red Blood Cells (RBC), which we are not interested in, and will not be analyzed further. Finally, we placed back some clusters that were found not to be in the correct CellLineage or CellType subset of cells. To do so, we followed 5 iterative steps, and in each directory we are providing one notebook detailing our choices made for each annotation step. We collected the annotations defined at the deepest resolution (64 populations, Level2 annotation) and group them into 15 major groups (Level 1 annotation, including Platelets and RBC).
  * *03_characterizing_CellTypes*: We run a differential Expression Analysis (DEA) to obtain the cell-type specific gene markers. We performed a comparison of our annotation with published PMBC atlases, to further validate our annotation strategy.
  * *SCGT00_CentralizedDataset*: Following the same structure as "Cell Annotation", but containing only scripts to define 15 major groups (Level 1 annotation, including Platelets and RBC) on the SCGT00 MAIN dataset.
  
3. **Downstream_Analysis**
In this folder we provide all the notebooks corresponding to the downstream analyses and the figures included in the manuscript.

    * *01_compositional_analaysis*: We applied [scCODA](https://github.com/theislab/scCODA) to assess difference in cell subpopulation among different inflammatory conditions. We considered *sex* and *age* as biological confounding factors, thus we corrected for them.
   
    * *02_gene_universe_definition*: This folder contains notebooks for defining the three gene sets i) Highly Variable, ii) Differentially Expressed, and iii) Manually Curated. Then, generate the scanpy objects needed for downstream analysis

    * *03_SEAcell_generation*: We applied [SEAcell](https://github.com/dpeerlab/SEACells) to aggregate cells into metacell. We applied this method to each patient indipendently, then we aggregated the results, using a [Nextflow](https://www.nextflow.io/) Pipeline. This step was needed in order to run Spectra.

    * *04_integration_with_annotation*: To obtain a batch-corrected embedding space we applied [scANVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html) considering cell annotation labels, as well as other methods such as [Harmony](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html), and [scGen](https://github.com/theislab/scgen). We, then sampled from the posterior learned during the training batch-corrected gene expression profiles, used in the following analyses. We compared results with non-corrected data, to be sure we are not biasing the main findings.

    * *05_SPECTRA*: We refined our curated list of immune-related genes, by employing [Spectra](https://github.com/dpeerlab/spectra), a tool able to identify a minimal set of genes related to specific functions in each cell-type (i.e., factors), starting from a predefined gene sets.

    * *06_inflammation_signatures*: This folder contains the necessary scripts used to compute and compare Spectra factors across diseases and cell types. First, we applied an signature scoring procedure (Univariate Linear Model (ULM) from [DecoupleR](https://decoupler-py.readthedocs.io/en/latest/generated/decoupler.run_ulm.html)). Then, we assess difference in signature activities between healthy and each inflammatory deseases, using a Linear Mixed Effect Model from [statsmodels](https://www.statsmodels.org/stable/mixed_linear.html) python library.

    * *07_gene_regulatory_network*: This directory includes the script to run GRN analyses that lead to identification of transcription factors associated with Interferon Response upregulation in SLE patients.

    * *08_gene_importance*: We applied a supervised classification approach, together with a post-hoc interpretability method, to allow the inference of the gene-wise importance, stratified by disease. We based our strategy on Gradient Boosted Decision Trees ([XGBoost](https://xgboost.readthedocs.io/en/stable/)). As GBDTs require post-hoc interpretability tools in order to infer explanations, we computed [SHAP](https://shap.readthedocs.io/en/latest/index.html) (SHapley Additive exPlanation) values. SHAP values explained the output of the classifier, in our case the predicted disease, as a sum of contributions of each feature, that is, the gene profiles. Such contributions correspond to the change in the expected model prediction when conditioning on the considered feature.

    * *09_patient_classifier*: We propose a novel computational pipeline  to exploit the cell embedding for classification of patients into inflammatory diseases, thus, turning the single-cell reference into a diagnostic tool. We first generated a cell type pseudobulk profile per patient by averaging the embedded features of the corresponding cells. Next, we trained independent classifier to assign correct disease labels, for each cell-type optimizing the hyperparameters with [Optuna](https://optuna.org/). To assess the performances of our approach, we proposed three scenarios i) a 5-fold cross validation strategy by splitting the full patient set into five balanced folds. ii) Prediction of unseen patients, randomly sampled from our CORE dataset and removed before starting the annotation task. iii) Prediction of unseen patients from new studies not included in our CORE dataset.
        * *SCGT00_CentralizedDataset*: Following the same structure as "09_patient_classifier", but containing only scripts referring to the patient classifier working with the SCGT00 dataset.
      

4. **Visualizing final Embedding Space**
This directory reports the notebooks needed for integrating the MAIN dataset and for generating the UMAP of the embedding space.
    * *SCGT00_CentralizedDataset*: Following the same structure as "Visualizing final Embedding Space", but containing only scripts referring to the final UMAP for the SCGT00 dataset.

6. **bin**
This directory includes several scripts with system requirements, enviornmental packages and versions, functions used in multiple scripts, some dictionaries used in the notebooks, as well as color palettes used for generating the figures.

6. **external_reference_data**
This directory includes external data required in the main analyses.


# Computational tools
The (most relevant) packages applied are listed below.
*Note: The requirements of many tools applied in this study were not compatible with each others. So, we generate different conda environments to execute specific analyses. We reported the list of libraries and tools with their respective versions at the end of each notebook.*

**Languages**
* Python > v.3.7
* R v.4.3.3

**Python tools:**
* [Decoupler-py](https://decoupler-py.readthedocs.io/en/latest/) v.1.6.0
* [Optuna](https://optuna.org/) v.3.6.0
* [scanpy](https://scanpy.readthedocs.io/en/stable/) v.1.9.8
* [scikit-learn](https://scikit-learn.org/stable/) v.1.4.1.post1
* [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/index_scrna.html) v.1.1.2
* [SEAcell](https://github.com/dpeerlab/SEACells) v.0.3.3
* [Spectra](https://github.com/dpeerlab/spectra) v.0.2.0
* [XGBoost](https://xgboost.readthedocs.io/en/stable/) py-xgboost-gpu: v2.0.3

**R tools:**
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) v.4.0.16
* [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) v.1.24.0

