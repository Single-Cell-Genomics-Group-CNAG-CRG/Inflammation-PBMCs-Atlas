
mainAdataFile = Channel.fromPath(params.mainAdataFilePath)

workflow {

    splitting_patients(mainAdataFile)
    splitting_patients.out
        .flatten()
        .map { [it.name - '.sample.h5ad', it] }
        .set { sampleID_h5ad }

    SEAcell_generation(sampleID_h5ad)

    merging_patient_SEAcell(mainAdataFile, SEAcell_generation.out.collect())

}

process splitting_patients {

    cpus 4
    memory '500 GB'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/inflammation_atlas_R1'

    storeDir params.outputdir + '/patient_split/'

    input:
        path mainAdataFile

    output:
        path "*.sample.h5ad"

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc

    # reading the main adata file

    adataM = sc.read_h5ad('${mainAdataFile}')
    
    if 'counts' in adataM.layers: #if false, we are assuming .X already contains the raw counts
        print('restoring counts')
        adataM.X = adataM.layers['counts'].copy()
        del adataM.layers

    del adataM.uns

    # removing unwanted cells

    excluded_CellType_str = "${params.exclude_CellTypes}"
    if excluded_CellType_str != '':
        excluded_CTs = excluded_CellType_str.split(',')
        annotationCol = "${params.annotation_col}"
        adataM = adataM[~adataM.obs[annotationCol].isin(excluded_CTs),:]
        print(f"{', '.join(excluded_CTs)} removed")

    # splitting by sampleID

    tot = adataM.obs.sampleID.unique().shape[0]

    for i, s_i in enumerate(adataM.obs.sampleID.unique()):
        
        adataS_i = adataM[adataM.obs.sampleID == s_i]

        print(f"Saving {s_i} ... ", end = '')
        adataS_i.write(f"{s_i}.sample.h5ad")
        print(f" done ({i+1}/{tot})")

    """
}

process SEAcell_generation {

    tag "${sampleID}"

    cpus 8

    time '11hours 59minutes'

    memory { task.attempt == 1 ? 16.GB : 64.GB }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 2

    storeDir params.outputdir + '/SEAcell_DFs/'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/SEAcell'

    input:
        tuple val(sampleID), path(sampleAdata)

    output:
        path "${sampleID}.seacell.pkl"


    """
    python ${projectDir}/SEAcellDefinition.py ${sampleAdata} ${params.SCxSEACell}
    """

}

process merging_patient_SEAcell {

    publishDir params.outputdir, mode: 'move', overwrite: true

    cpus 4
    memory '64 GB'

    conda '/scratch_isilon/groups/singlecell/shared/conda_env/inflammation_atlas_R1'

    input:
        path mainAdata
        path '*.seacell.pkl'

    output:
        path params.outputFileName

    script:
    """
    #!/usr/bin/env python
    from glob import glob
    import scanpy as sc
    import pandas as pd

    # reading main adata file
    # adataM = sc.read_h5ad('${mainAdata}')

    # collecting sample adata files
    samplePklFileL = glob('*.seacell.pkl')
    print(f"Collected {len(samplePklFileL)} sample dataframes")

    # read and merge them
    samplePklList = []
    for s_i in samplePklFileL:
        samplePklList.append(pd.read_pickle(s_i))

    seacellDF = pd.concat(samplePklList, axis = 0, ignore_index = False, verify_integrity = True)

    seacellDF.to_pickle('${params.outputFileName}')

    # adding SEAcells aggregation to the main adata file
    # adataM_obs = adataM.obs.merge(seacellDF, left_index = True, right_index = True, how = 'left')

    
    """

}