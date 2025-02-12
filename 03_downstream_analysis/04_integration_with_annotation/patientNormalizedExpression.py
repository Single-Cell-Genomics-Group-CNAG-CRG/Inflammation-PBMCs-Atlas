import sys
import os

import scanpy as sc
import pandas as pd

import decoupler

import torch
torch.set_float32_matmul_precision('high')

# Set random seed
random_seed = 5

import scvi
scvi.settings.dl_num_workers = 8
scvi.settings.seed = random_seed

import torch
torch.set_float32_matmul_precision('high')
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
torch.multiprocessing.set_sharing_strategy('file_system')


mainAdataFilePath = sys.argv[1]
scANVImodelDirPath = sys.argv[2]
SEAcellDFpath = sys.argv[3]

sampleID = sys.argv[4]
annotationCol = sys.argv[5]
exclude_CellTypes = sys.argv[6].split(',')
batches_list = sys.argv[7].split(',')

print(sampleID)


os.mkdir('./patient_seacell')
os.mkdir('./patient_pseudobulk')
os.mkdir('./patient_cellType')



def aggregate_cells(adata = None,
                    SEAcellDF = None,
                    feature_mtx = None, 
                    sparse = True,
                    SEAcellIDcol = 'SEACell', 
                    how = 'sum', 
                    return_adata = True,
                    annotationCol = None, 
                    returnAllCellTypeProportion = False,
                    SEACellMetadataCols = [],
                    workers = 16):

################ function ################################
    def process_cell(SC_i):
        import numpy as np
        cellIdx_i = adata_obs['SEACell'] == SC_i
        if howL in ['mean','average']:
            aggr_array_i = np.asarray(np.mean(data[cellIdx_i,:], axis = 0)).reshape(-1)
        elif howL == 'sum':
            aggr_array_i = np.asarray(np.sum(data[cellIdx_i,:], axis = 0)).reshape(-1)
        else:
            aggr_array_i = []
            
        if computeCTfreq:
            # compute the most abundat cell type and the related proportion
            CTproportion_i = annotations.loc[cellIdx_i].value_counts(normalize=True, sort = True, ascending = False)
            MostAbCT_i = CTproportion_i.index[0][0]
            MostAbProp_i = CTproportion_i.iloc[0]
            nSC_i = cellIdx_i.sum()
            if returnAllCTProp:
                return(SC_i, aggr_array_i, MostAbCT_i, MostAbProp_i, nSC_i, CTproportion_i)
            return (SC_i, aggr_array_i, MostAbCT_i, MostAbProp_i, nSC_i)
        return (SC_i, aggr_array_i)
###################################################3
    
    assert ((adata is None) != (SEAcellDF is None))

    if adata is None:
        how = None
        return_adata = False
        
    
    import scanpy as sc
    from scipy.sparse import csr_matrix
    import pandas as pd
    import numpy as np
    from tqdm import tqdm

    if workers > 1:
        from multiprocess import Pool

    global data
    global adata_obs
    global computeCTfreq
    global returnAllCTProp

    computeCTfreq = False
    returnAllCTProp = False
    
    if annotationCol is not None:
        global annotations
        if adata is None:
            annotations = SEAcellDF[[annotationCol]]
        else:
            annotations = adata.obs[[annotationCol]]
        computeCTfreq = True
        adata_obs = SEAcellDF[[SEAcellIDcol] + SEACellMetadataCols + [annotationCol]]
    else:
        adata_obs = adata.obs[[SEAcellIDcol] + SEACellMetadataCols]

    if how is not None:
        if how not in ['sum','mean','average']:
            raise ValueError(f" {how} not supported. Choose one between 'sum' or 'mean'")        
        
    global howL 
    if how is None:
        howL = None
    else:
        howL = how.lower()
        
    if how is not None:
        if (feature_mtx.lower() == 'x') or (feature_mtx is None):
            data = adata.X
            featureName = adata.var_names
        elif feature_mtx.lower() == 'raw':
            if adata.raw is None:
                raise ValueError("RAW data are not present")
            else:
                data = adata.raw.X
                featureName = adata.raw.var_names
        elif feature_mtx in adata.layers:
            data = adata.layers[feature_mtx]
            featureName = adata.var_names
        elif feature_mtx in adata.obsm:
            data = adata.obsm[feature_mtx]
            featureName = None
            if return_adata:
                print("WARNING: var metadata could not correspond to features aggregated.")
                print("If needed, specify return_adata = False to return only the aggregated matrix.")
        else:
            raise ValueError(f" {feature_mtx} matrix not found")
    

    SEACellID = adata_obs[SEAcellIDcol].unique().astype(str)

    aggregationList = []
    if annotationCol is not None:
        MostAbCT_list = []
        MostAbProp_list = []
        nSC_list = []
        if returnAllCellTypeProportion:
            returnAllCTProp = True
            CTproDF_list = []            

    if workers == 1:
        for SC_i in tqdm(SEACellID):
            cellIdx_i = adata_obs['SEACell'] == SC_i
            if how in ['mean','average']:
                aggregationList.append(np.asarray(np.mean(data[cellIdx_i,:], axis = 0)).reshape(-1))
            elif how == 'sum':
                aggregationList.append(np.asarray(np.sum(data[cellIdx_i,:], axis = 0)).reshape(-1))
                
            if annotationCol is not None:
                # compute the most abundat cell type and the related proportion
                CTproportion_i = adata_obs.loc[cellIdx_i, annotationCol].value_counts(normalize=True, sort = True, ascending = False)
                MostAbCT_list.append(CTproportion_i.index[0])
                MostAbProp_list.append(CTproportion_i.iloc[0])
                nSC_list.append(cellIdx_i.sum())
                if returnAllCTProp:
                    CTproDF_list.append(CTproportion_i)
    else:
        p = Pool(workers)
        resultList = list(tqdm(p.imap(process_cell, SEACellID), total=len(SEACellID)))
        p.close()
        SEACellID = [] #be sure that the results are in the correct order 
        print('Collecting results from workers')
        for t in tqdm(resultList):
            SEACellID.append(t[0])
            aggregationList.append(t[1])
            if computeCTfreq:
                MostAbCT_list.append(t[2])
                MostAbProp_list.append(t[3])
                nSC_list.append(t[4]) 
                if returnAllCTProp:
                    CTproDF_list.append(t[5])
    if how is not None:
        SEAcellGroupedDF = np.vstack(aggregationList)
    else:
        SEAcellGroupedDF = None

    if returnAllCTProp:
        CTpropDF = pd.concat(CTproDF_list, axis=1).reset_index().set_index(annotationCol)
        CTpropDF.columns = SEACellID
        CTpropDF = CTpropDF.transpose()
        CTpropDF.columns = CTpropDF.columns.tolist()
        
        
    if len(SEACellMetadataCols) > 0:
        SEACellMD = adata_obs[[SEAcellIDcol] + SEACellMetadataCols].drop_duplicates().set_index(SEAcellIDcol)

    # building obs dataframe
            # storing SEACell IDs
    SEACellAdata_obs = pd.DataFrame(index=SEACellID)
    if annotationCol is not None:
        SEACellAdata_obs['MostAbundantCellType'] = MostAbCT_list
        SEACellAdata_obs['MostAbundantCellType_proportion'] = MostAbProp_list
        SEACellAdata_obs['n_single_cell'] = nSC_list

    if len(SEACellMetadataCols) > 0:
        SEACellAdata_obs = SEACellAdata_obs.merge(SEACellMD, left_index = True, right_index = True, how = 'left')

    if return_adata:
        # creating the SEACell aggregated anndata object
        if sparse:
            SEACellAdata = sc.AnnData(csr_matrix(SEAcellGroupedDF), 
                                      obs = SEACellAdata_obs,
                                      var = adata.var.copy(deep=True))#, dtype=csr_matrix(SEAcellGroupedDF).dtype)
        else:
            SEACellAdata = sc.AnnData(SEAcellGroupedDF, 
                                      obs = SEACellAdata_obs,
                                      var = adata.var.copy(deep=True))#, dtype=csr_matrix(SEAcellGroupedDF).dtype)            

        if returnAllCTProp:
            SEACellAdata.obsm['CellType_proportions'] = CTpropDF

        return SEACellAdata

    else:
        # returning the aggregated matrix with same col namesas the original one
        if how is not None:
            SEAcellGroupedDF = pd.DataFrame(SEAcellGroupedDF)
            if featureName is not None:
                SEAcellGroupedDF.columns = featureName
        
            # storing SEACell IDs as row index 
            SEAcellGroupedDF.set_index(SEACellID, inplace=True)
        
            if returnAllCTProp:
                return SEAcellGroupedDF.sort_index(), SEACellAdata_obs.sort_index(), CTpropDF.sort_index()
            else:
                return SEAcellGroupedDF.sort_index(), SEACellAdata_obs.sort_index()
        else:
            if returnAllCTProp:
                return SEACellAdata_obs.sort_index(), CTpropDF.sort_index()
            else:
                return SEACellAdata_obs.sort_index()
            
#### Starting of the script

# Load the h5ad file
adata = sc.read_h5ad(mainAdataFilePath)# ,backed='r+', chunk_size=50000)
adata.obs['binned_age'] = adata.obs['binned_age'].astype(str)


#loading scANVI model
scanvi_model = scvi.model.SCANVI.load(scANVImodelDirPath, adata=adata) 

#loading SEAcell dataframe
SEAcellDF = pd.read_pickle(SEAcellDFpath)


#GET NORMALIZED EXPRESSION
adata_obs = adata.obs.reset_index().copy()
adata_obsFilt = adata_obs.loc[(~adata_obs[annotationCol].isin(exclude_CellTypes)) & (adata_obs['sampleID'] == sampleID)]
indices = adata_obsFilt.index.tolist()

normalized_expression_S_CT = scanvi_model.get_normalized_expression(n_samples=25, 
                                                                    transform_batch = ['3_GEX_V3','5_GEX_V2','3_GEX_V2','5_GEX_V1'],
                                                                    indices = indices,
                                                                    library_size=10000
                                                                   )

# creating adata from normalized expression
sampleNormAdata = sc.AnnData(X = normalized_expression_S_CT, 
                             obs = adata_obsFilt.set_index('cellID'))


# generating pseudobulk
print()
pseudobulk_i = decoupler.get_pseudobulk(adata=sampleNormAdata, 
                                        min_cells=0, 
                                        sample_col = 'sampleID',
                                        groups_col=annotationCol, 
                                        layer=None,#'counts', 
                                        mode='mean')

pseudobulk_i.obs.drop(['_scvi_batch','_scvi_labels'], axis=1, inplace=True)

# remove layers (it includes only the proportion of cells that express each gene. 
del pseudobulk_i.layers

pseudobulk_i.write_h5ad(f"./patient_pseudobulk/{sampleID}.pseudobulk.h5ad", compression='gzip')

# Cycle over cell-types to save a splitted anndata
sampleNormAdata.obs.drop(['_scvi_batch', '_scvi_labels'], axis=1, inplace=True)

for ct in sampleNormAdata.obs[annotationCol].unique():
    print(ct)
    ad_i = sampleNormAdata[sampleNormAdata.obs[annotationCol] == ct].copy()
    # sc.pp.log1p(ad_i)
    #ad_i.write_h5ad(f"{sampleID}_cellType_{ct}.log1p.h5ad", compression='gzip')
    ad_i.write_h5ad(f"./patient_cellType/{sampleID}_cellType_{ct}.h5ad", compression='gzip')

### Aggregate SEAcell
# We are computing the average of the normalized data. Then, we log1p scale them

sampleNormAdata.obs = sampleNormAdata.obs.drop(annotationCol, axis=1).merge(SEAcellDF.drop(['disease','sampleID'], axis=1), right_index=True, left_index=True, how='left')

# Keeping only cells that were included in one SEAcell
sampleNormAdataFilt = sampleNormAdata[~sampleNormAdata.obs.SEACell.isna()]

if sampleNormAdataFilt.shape[0] > 0:

  SEAcellNorm = aggregate_cells(adata = sampleNormAdataFilt,
                      SEAcellDF = None,
                      feature_mtx = 'X', 
                                sparse = False,
                      SEAcellIDcol = 'SEACell', 
                      how = 'mean', 
                      return_adata = True,
                      annotationCol = None, 
                      returnAllCellTypeProportion = False,
                      SEACellMetadataCols = [annotationCol,'sampleID','disease'],
                      workers = 1)

  #SEAcellNorm.write_h5ad(f"{sampleID}.seacell.log1p.h5ad", compression='gzip')
  SEAcellNorm.write_h5ad(f"./patient_seacell/{sampleID}.seacell.h5ad", compression='gzip')