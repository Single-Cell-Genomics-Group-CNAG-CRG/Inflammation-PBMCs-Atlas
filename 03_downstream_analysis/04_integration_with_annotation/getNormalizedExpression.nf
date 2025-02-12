
params.mainAdataFilePath = '/scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/02_gene_universe_definition/results/04_MAIN_geneUniverse_noRBCnPlatelets.h5ad'

params.SEAcellDFpath = '/scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/03_SEAcell_generation/results/01_MAIN_SEAcell_annotated_LowPurityFilt.pkl'

params.scANVImodelDirPath = '/scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/04_integration_with_annotation/results/scANVI_model_fineTuned_lowLR_noRBCnPlat/'

params.outputdir = '/scratch_isilon/groups/singlecell/shared/projects/Inflammation-PBMCs-Atlas/03_downstream_analysis/04_integration_with_annotation/results/normalized_adatas_nextflow_v3'

params.exclude_CellTypes = 'RBC,Platelets'
params.annotation_col = 'Level1'
params.sample_batches = '3_GEX_V2' 

mainAdataFile = Channel.fromPath(params.mainAdataFilePath)

SEAcellDFfile = Channel.fromPath(params.SEAcellDFpath)

scANVImodelFile = Channel.fromPath(params.scANVImodelDirPath)

workflow {

    getSampleList(mainAdataFile)
    getSampleList.out
        .splitCsv( header: true )
        .map { row -> "${row.sampleID}" }
        .set { sampleList }

    getNormExp_aggrSEAcellSample_splitCellType(sampleList, mainAdataFile.first(), scANVImodelFile.first(), SEAcellDFfile.first())
    getNormExp_aggrSEAcellSample_splitCellType.out.celltype_adatas
            .flatten()
            .map {it ->tuple(it.toString().split("_cellType_")[1].split(".h5ad")[0], it)}
            .groupTuple()
            .set { sample_cellType_file }

    mergeCellType(sample_cellType_file)
    mergePseudobulk(getNormExp_aggrSEAcellSample_splitCellType.out.pseudobulk_adatas.collect())
    mergeSEAcell(getNormExp_aggrSEAcellSample_splitCellType.out.seacell_adatas.collect())
}

process getSampleList {

    cpus 4
    memory '32 GB'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/inflammation_atlas_R1'

    storeDir params.outputdir

    input:
        path mainAdataFile

    output:
        path "sampleList.csv"


    """
    #!/usr/bin/env python

    import scanpy as sc
    import pandas as pd

    # reading the main adata file

    adataM = sc.read_h5ad('${mainAdataFile}')
    
    excluded_CellType_str = "${params.exclude_CellTypes}"
    if excluded_CellType_str != '':
        excluded_CTs = excluded_CellType_str.split(',')
        adataM = adataM[~adataM.obs["${params.annotation_col}"].isin(excluded_CTs),:]
        print(f"{', '.join(excluded_CTs)} removed")   

    sampleInfo = pd.DataFrame(adataM.obs.groupby('sampleID', observed=True).size())
    sampleInfo.columns = ['n_cells']
    sampleInfo.to_csv('./sampleList.csv')
    """
}

process getNormExp_aggrSEAcellSample_splitCellType {

    tag "${sampleID}" 

    memory { task.attempt == 1 ? 64.GB : 128.GB }
    cpus { task.attempt == 1 ? 16 : 32 }

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/scvi-v112'

    storeDir params.outputdir

    // publishDir params.outputdir + '/sample_seacell', saveAs: { (it ==~ /.*.seacell.log1p.h5ad/) ? "$it" : null }, mode: 'copy'
    // publishDir params.outputdir + '/sample_cellType', saveAs: { (it ==~ /.*_cellType_.*/) ? "$it" : null }, mode: 'copy'

    input:
        val sampleID
        path mainAdataFile
        path scANVImodelFile
        path SEAcellDFfile

    output:
        path "patient_seacell/${sampleID}.seacell.h5ad", emit: seacell_adatas//, optional: false
        path "patient_pseudobulk/${sampleID}.pseudobulk.h5ad", emit: pseudobulk_adatas//, optional: false
        path "patient_cellType/${sampleID}_cellType_*.h5ad", emit: celltype_adatas  


    """
    python ${projectDir}/patientNormalizedExpression.py '${mainAdataFile}' \
                                                        '${scANVImodelFile}' \
                                                        '${SEAcellDFfile}' \
                                                        '${sampleID}' \
                                                        '${params.annotation_col}' \
                                                        '${params.exclude_CellTypes}' \
                                                        '${params.sample_batches}'
    """
}

process mergeCellType {

    tag "${cellType}" 

    cpus 16
    memory '500 GB'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/inflammation_atlas_R1'

    storeDir params.outputdir + '/cellType_adata_merged'

    input:
        tuple val(cellType), path("cellTypeAdatas/${cellType}_*.h5ad")

    output:
        path "${cellType}_adataMerged.log1p.h5ad"

    """
    #!/usr/bin/env python

    from glob import glob
    import scanpy as sc
    import pandas as pd
    import numpy as np

    # collecting sample adata files
    adatas_files = glob('cellTypeAdatas/${cellType}_*.h5ad')
    print(f"Collected {len(adatas_files)} sample seacell adatas")

    # read and merge them
    adata_list = []
    for s_i in adatas_files:
        adata_list.append(sc.read_h5ad(s_i))

    adataMerged = sc.concat(adata_list, join='outer', fill_value=-1, index_unique=None)

    # log-scaling
    sc.pp.log1p(adataMerged)

    print(np.min(adataMerged.X))

    adataMerged.write_h5ad('${cellType}_adataMerged.log1p.h5ad', compression='gzip')

    """

}

process mergePseudobulk {
    cpus 16
    memory '100 GB'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/inflammation_atlas_R1'

    storeDir params.outputdir

    input:
        path 'pseudobulkAdatas/sampleID_*.h5ad'

    output:
        path "pseudobulkAdataMerged.h5ad"

    """
    #!/usr/bin/env python

    from glob import glob
    import scanpy as sc
    import pandas as pd
    import numpy as np

    # collecting sample adata files
    pseudobulk_adatas_files = glob('pseudobulkAdatas/sampleID_*.h5ad')
    print(f"Collected {len(pseudobulk_adatas_files)} sample adatas")

    # read and merge them
    pseudobulk_adata_list = []
    for s_i in pseudobulk_adatas_files:
        pseudobulk_adata_list.append(sc.read_h5ad(s_i))

    pseudobulkAdataMerged = sc.concat(pseudobulk_adata_list, join='outer', fill_value=-1, index_unique=None)

    pseudobulkAdataMerged.write_h5ad('pseudobulkAdataMerged.h5ad', compression='gzip')

    """
}


process mergeSEAcell {

    cpus 16
    memory '500 GB'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/inflammation_atlas_R1'

    storeDir params.outputdir

    input:
        path 'SEAcellAdatas/sampleID_*.h5ad'

    output:
        path "SEAcellAdataMerged.log1p.h5ad"

    """
    #!/usr/bin/env python

    from glob import glob
    import scanpy as sc
    import pandas as pd
    import numpy as np

    # collecting sample adata files
    seacell_adatas_files = glob('SEAcellAdatas/sampleID_*.h5ad')
    print(f"Collected {len(seacell_adatas_files)} sample seacell adatas")

    # read and merge them
    seacell_adata_list = []
    for s_i in seacell_adatas_files:
        seacell_adata_list.append(sc.read_h5ad(s_i))

    seacellAdataMerged = sc.concat(seacell_adata_list, join='outer', fill_value=-1, index_unique=None)

    sc.pp.log1p(seacellAdataMerged)

    seacellAdataMerged.write_h5ad('SEAcellAdataMerged.log1p.h5ad', compression='gzip')

    """
}