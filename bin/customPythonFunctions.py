def getClassMetricsDF(data=None, disease_list = None, y_true_cn = 'disease', y_pred_cn = 'disease_pred', cell_type_cn = 'Level1', include_aggregated = True):

    import pandas as pd
    from sklearn.metrics import recall_score, precision_score, f1_score, accuracy_score

    classMtx = dict()
    for f in [recall_score, precision_score, f1_score, accuracy_score]:
        classMtx[f.__name__] = dict()
        for ct in data[cell_type_cn].unique():
            if f.__name__ == 'accuracy_score':
                disease_res = []
                for d in disease_list:
                    disease_res.append(
                        f(y_true = data.query(f"{cell_type_cn} == @ct")[y_true_cn] == d, 
                          y_pred = data.query(f"{cell_type_cn} == @ct")[y_pred_cn] == d
                         )                        
                    )
                    
                classMtx[f.__name__][ct] = disease_res
            else:
                classMtx[f.__name__][ct] = f(y_true = data.query(f"{cell_type_cn} == @ct")[y_true_cn],
                                             y_pred = data.query(f"{cell_type_cn} == @ct")[y_pred_cn], 
                                             labels = disease_list, 
                                             average = None)

        if include_aggregated:
            # computing the metric on the overall result, not stratified by cell-types
            if f.__name__ == 'accuracy_score':
                disease_res = []
                for d in disease_list:
                    disease_res.append(
                        f(y_true = data[y_true_cn] == d, 
                          y_pred = data[y_pred_cn] == d)                        
                    )
                classMtx[f.__name__]['aggregated'] = disease_res
            else:
                classMtx[f.__name__]['aggregated'] = f(y_true = data[y_true_cn], 
                                                       y_pred = data[y_pred_cn], 
                                                       labels = disease_list, average = None)

        df = pd.concat([pd.DataFrame.from_dict(classMtx[f]).T.assign(metric=f) for f in classMtx]).reset_index()
        df.columns = ['stratification'] + disease_list + ['metric']
        
    return df.loc[:,['metric', 'stratification'] + disease_list]


def generate_shap_data(
    cell_type: str = '',
    shap_stats_path: str = '', 
    adata_path: str = '',
    gene_symbol_df_path: str='',
    stat: str = 'mean_abs',
    category_col: str = 'disease',
    expressed_gene_cellTypes_path: str = ''
):
    
    import numpy as np
    import anndata as ad
    import pandas as pd
    
    shap_stats = np.load(shap_stats_path)[stat]

    symbol_df = pd.read_pickle(gene_symbol_df_path)

    adata = ad.read_h5ad(adata_path, backed = 'r')

    adata.var = adata.var.merge(symbol_df, left_index=True, right_index=True, how='left')

    genes = adata.var['symbol'].values

    categories = adata.obs[category_col].cat.categories

    ## pandas dataframe

    shapDF = pd.DataFrame(shap_stats, index=genes, columns=categories)

    ## keeping only selected genes
    well_expressed_symbols = pd.read_csv(expressed_gene_cellTypes_path).query("CellType == @cell_type and `% cells` > 5").symbol.values

    shapDF_filt = shapDF.iloc[shapDF.index.isin(well_expressed_symbols),:]

    return shapDF_filt
    
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
            CTproportion_i = annotations[cellIdx_i].value_counts(normalize=True, sort = True, ascending = False)
            MostAbCT_i = CTproportion_i.index[0]
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
            annotations = SEAcellDF[annotationCol].cat.remove_unused_categories()
        else:
            annotations = adata.obs[annotationCol].cat.remove_unused_categories()
        computeCTfreq = True
        if adata is None:
            adata_obs = SEAcellDF[[SEAcellIDcol] + SEACellMetadataCols + [annotationCol]]
        else:
            adata_obs = adata.obs[[SEAcellIDcol] + SEACellMetadataCols + [annotationCol]]
    else:
        if adata is None:
            adata_obs = SEAcellDF[[SEAcellIDcol] + SEACellMetadataCols]
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
                CTproportion_i = annotations[cellIdx_i].value_counts(normalize=True, sort = True, ascending = False)
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

def plot_clusters_separately(data, leiden_col, title='', basis="umap",size=7):
    
    clusters = data.obs[leiden_col].cat.categories
    n_clust = len(clusters)
    row_num = np.ceil(n_clust/4).astype(int)
    fig, ax = plt.subplots(nrows=row_num, ncols=4, figsize=(10, row_num*3))
    ax = ax.ravel()

    plt.rcParams.update({'axes.labelsize': 0})
    for clust in range(0, len(clusters)):
        # print(clust)
        clust_name = clusters[clust]
        sc.set_figure_params(vector_friendly=True, dpi_save=300)  # Makes PDFs of scatter plots much smaller in size but still high-quality
        sc.pl.scatter(data, 
                      color=leiden_col, 
                      basis=basis, 
                      palette=["red"], 
                      groups=str(clust_name), 
                      size=size, 
                      show=False, 
                      ax=ax[clust], 
                      legend_loc="none", 
                      title=leiden_col)
        ax[clust].set_title("Cluster " + clust_name, size=15)
        ax[clust].xaxis.label.set_size(0)
        ax[clust].yaxis.label.set_size(0)
        ax[clust].set_rasterization_zorder(10)

    # delete empty axes
    sub_plots_no = row_num*4
    if (n_clust < sub_plots_no):
        for i in range(1,sub_plots_no-n_clust+1):
            fig.delaxes(ax.flatten()[sub_plots_no-i])
    # if (pdf==None):
    #     pdf.savefig(fig)
    # else:
    #     #plt.show()
    #     pdf.savefig(fig, bbox_inches="tight")
    #     #plt.savefig(pdf)
        
    fig.suptitle(title)
    plt.show()

def plot_UMAPs_lineages(adata, marker_genes=None, path_figure = None, saveFig=False, dpi_fig_save=None):
    #global adata
    import scanpy as sc
    import matplotlib.pyplot as plt

    random_indices = balanced_sample(adata.obs, cols = ["sampleID"], frac = 0.1, shuffle = True, random_state = 42).cellID
    nrow=len(marker_genes)
    ncol=max([len(vs) for vs in marker_genes.values()])
    fig,axs=plt.subplots(nrow,ncol,figsize=(3*ncol,3*nrow))
    
    # Plot expression for every marker on the corresponding Axes object
    for row_idx,(cell_type,markers) in enumerate(marker_genes.items()):
        col_idx=0
        for marker in markers:
            if marker not in adata.var_names:
                continue
            ID2SymbolDF = generateID2SymbolDF(varDF = adata.var, 
                                              symbolList = marker, 
                                              ID_col = 'index', symbols_col = 'symbol', 
                                              HUGOstatus_col = 'HUGO_status', 
                                              behaviour = 'all')
            ax=axs[row_idx,col_idx]
            sc.pl.embedding(basis = 'X_umap_scVI', 
                            adata = adata[random_indices, :], 
                            color=ID2SymbolDF["gene_id"], title= ID2SymbolDF["symbol"],
                            ax = ax,show=False,
                            frameon=False,
                            vmin="p1", vmax="p99",
                            s= 0.15, 
                            use_raw = False)
            # Add cell type as row label - here we simply add it as ylabel of
            # the first Axes object in the row
            if col_idx==0:
                # We disabled axis drawing in UMAP to have plots without background and border
                # so we need to re-enable axis to plot the ylabel
                ax.axis('on')
                ax.tick_params(
                    top='off', bottom='off', left='off', right='off',
                    labelleft='on', labelbottom='off')
                ax.set_ylabel(cell_type+'\n', rotation=90, fontsize=14)
                ax.set_xlabel('')
                ax.set(frame_on=False)
            col_idx+=1
        # Remove unused column Axes in the current row
        while col_idx<ncol:
            axs[row_idx,col_idx].remove()
            col_idx+=1
            
    if (saveFig == True) and ~(path_figure is None) :
        plt.savefig(path_figure, bbox_inches='tight', pad_inches=0, dpi=dpi_fig_save)
        
    return None
    
def generateID2SymbolDF(varDF = None, symbolList = None, ID_col = 'index', symbols_col = 'symbol', HUGOstatus_col = 'HUGO_status', behaviour = 'all'):

    import pandas as pd
    # behaviour defines what happen if there are duplicated symbol 
    # 'all': returned list will include all the corresponding IDs
    # 'official': prioritize gene ID with official status. If none of the symbols are official all are returned  

    gene_namesDF = varDF.loc[varDF[symbols_col].isin(symbolList)]

    gene_namesDF[symbols_col] = pd.Categorical( gene_namesDF[symbols_col], categories=symbolList, ordered=True)

    if gene_namesDF[symbols_col].duplicated().sum() == 0:

        if ID_col == 'index':
            return gene_namesDF.reset_index(names = 'gene_id')[['gene_id',symbols_col]].sort_values(symbols_col)
        else:
            return gene_namesDF[[ID_col,symbols_col]].sort_values(symbols_col)
    else:

        duplicatedGeneList = gene_namesDF[symbols_col][gene_namesDF[symbols_col].duplicated()].tolist()
        
        print("The following gene symbols are duplicated: ")
        [print(f"{g},") for g in duplicatedGeneList]

        if behaviour == 'all':
            if ID_col == 'index':
                return gene_namesDF.reset_index(names = 'gene_id')[['gene_id',symbols_col]].sort_values(symbols_col)
            else:
                return gene_namesDF[[ID_col,symbols_col]].sort_values(symbols_col)
                
        elif behaviour == 'official':
            if HUGOstatus_col is None:
                raise ValueError("Please, specify column with HUGO status")
            else:
                uniqueSymbolDF = gene_namesDF.loc[~gene_namesDF[symbols_col].duplicated(keep=False)]
                duplicatedSymbolDF = gene_namesDF.loc[gene_namesDF[symbols_col].duplicated(keep=False)]

            officialDF = duplicatedSymbolDF.loc[duplicatedSymbolDF[HUGOstatus_col] == 'official']
            nonOfficial = duplicatedSymbolDF.loc[~duplicatedSymbolDF[symbols_col].isin(officialDF[symbols_col])]

            gene_namesDF = pd.concat([uniqueSymbolDF,officialDF,nonOfficial], axis = 0)
            
            if ID_col == 'index':
                return gene_namesDF.reset_index(names = 'gene_id')[['gene_id',symbols_col]].sort_values(symbols_col)
            else:
                return gene_namesDF[[ID_col,symbols_col]].sort_values(symbols_col)
    


def balanced_sample(df, cols = None, n = None, frac = None, shuffle = False, random_state = 42):

    import pandas as pd

    if not ((n is None) != (frac is None)): 
        print("Error: please specify n or frac, not both")
        return None
        
    # Group by the columns and apply the sample function

    if cols is None:
        df_sampled = df.apply(lambda x: x.sample(n = n, frac=frac, replace = False, random_state=random_state))
    else:
        df_sampled = df.groupby(cols, observed=True).apply(lambda x: x.sample(n = n, frac=frac, replace = False, random_state=random_state))
        df_sampled = df_sampled.drop(cols, axis=1, errors='ignore').reset_index()

    if shuffle:
        return df_sampled.sample(frac=1, random_state=random_state)
    else:
        return df_sampled
 
def get_deltas_v2(scGenModel, batchLevels = []):
    from scvi import REGISTRY_KEYS
    from anndata import AnnData
    import pandas as pd
    from tqdm import tqdm
    from collections import defaultdict
    import numpy as np
    
    """
    Modified function from scGen Removes batch effects.
    https://scgen.readthedocs.io/en/stable/_modules/scgen/_scgen.html#SCGEN.batch_removal
    """
    latent_all = scGenModel.get_latent_representation(adata=None)
    
    adata_latent = AnnData(latent_all)
    adata_latent.obs = scGenModel.adata.obs.copy(deep=True)

    batchDF = pd.DataFrame(adata_latent.obs.groupby(batchLevels).size())
    batchDF.columns = ['count']
    batchDF = batchDF[batchDF['count'] > 0]
    batchDF.reset_index(inplace=True)
    
    if len(batchLevels) < 2:
        maxBatches = batchDF.sort_values('count', ascending=False).iloc[0,0]
        batchDF[f"{batchLevels[0]}_max"] = maxBatches
    else:
        maxBatches = batchDF[batchDF.groupby(batchLevels[0:-1])['count'].transform(max) == batchDF['count']]
        maxBatches.columns = ['{}_max'.format(c) for c in maxBatches.columns]
        maxBatches.drop('count_max', axis = 1, inplace=True)
        batchDF = batchDF.merge(maxBatches, how='left', left_on=batchLevels[0:-1], right_on=['{}_max'.format(c) for c in batchLevels[0:-1]])    

    nested_dict = lambda: defaultdict(nested_dict)
    deltas = nested_dict()

    for i in tqdm(range(batchDF.shape[0])):#[batchLevels].iterrows():
        Bcomb = batchDF.iloc[i,:]
        subsetObs = adata_latent.obs[batchLevels]
        subsetObsMax = adata_latent.obs[batchLevels]
    
        tempdeltas = deltas
        
        for b_ in batchLevels:
            v_ = Bcomb[b_]
            q_ = f"{b_} == '{v_}'"  
            subsetObs = subsetObs.query(q_)
            
            v_max = Bcomb['{}_max'.format(b_)]
            q_max = f"{b_} == '{v_max}'" 
            subsetObsMax = subsetObsMax.query(q_max)
            
            if b_ == batchLevels[-1]:
                # compute delta
                delta = np.average(adata_latent[adata_latent.obs.index.isin(subsetObsMax.index)].X, axis=0) - np.average(adata_latent[adata_latent.obs.index.isin(subsetObs.index)].X, axis=0)
                tempdeltas[v_] = delta
            else:
                tempdeltas = tempdeltas[v_]

    return deltas


def batch_removal_v2(scGenModel, adata=None, deltas=None, batchLevels = []):
    """
    Modified function from scGen Removes batch effects.
    https://scgen.readthedocs.io/en/stable/_modules/scgen/_scgen.html#SCGEN.batch_removal
    """
    from anndata import AnnData
    import pandas as pd
    from tqdm import tqdm
    from collections import defaultdict
    import torch
    import numpy as np
    
    latent_all = scGenModel.get_latent_representation(adata=adata)

    
    adata_latent = AnnData(latent_all)
    
    if adata is None:
        adata = scGenModel.adata
    
    adata_latent.obs = adata.obs.copy(deep=True)

    batchDF = pd.DataFrame(adata_latent.obs.groupby(batchLevels).size())
    batchDF.columns = ['count']
    batchDF = batchDF[batchDF['count'] > 0]
    batchDF.reset_index(inplace=True)

    for i in tqdm(range(batchDF.shape[0])):#[batchLevels].iterrows():
        Bcomb = batchDF.iloc[i,:]
        subsetObs = adata_latent.obs[batchLevels]
        
        tempdeltas = deltas
        
        for b_ in batchLevels:
            v_ = Bcomb[b_]
            q_ = f"{b_} == '{v_}'"  
            subsetObs = subsetObs.query(q_)
            
            if v_ not in tempdeltas:
                continue
            
            if isinstance(tempdeltas[v_], np.ndarray):
                # compute delta
                adata_latent[adata_latent.obs.index.isin(subsetObs.index)].X  += tempdeltas[v_] #+ adata_latent[adata_latent.obs.index.isin(subsetObs.index)].X 
            else:
                tempdeltas = tempdeltas[v_]
    with torch.no_grad():                
        corrected = AnnData(
            scGenModel.module.generative(torch.Tensor(adata_latent.X))["px"]
            .cpu()
            .numpy(),
            obs=adata_latent.obs,
        )
    corrected.var_names = adata.var_names.tolist()
    corrected = corrected[adata.obs_names]
    if adata.raw is not None:
        adata_raw = AnnData(X=adata.raw.X, var=adata.raw.var)
        adata_raw.obs_names = adata.obs_names
        corrected.raw = adata_raw
    corrected.obsm["latent"] = adata_latent.X
    corrected.obsm["corrected_latent"] = scGenModel.get_latent_representation(
        corrected
    )
    return corrected
    
def cumulative_explained_variance(adata, expVarThr=None, n_pc=None, figsize = (8,8)):
    import numpy as np
    import matplotlib.pyplot as plt

    cs = np.cumsum(adata.uns['pca']['variance'] / np.sum(adata.uns['pca']['variance']))

    if figsize is not None:
        plt.figure(figsize=figsize)
        plt.plot(range(1,adata.uns['pca']['variance'].shape[0]+1), 
                    np.cumsum(adata.uns['pca']['variance'] / np.sum(adata.uns['pca']['variance']))
                   )

    res = dict({'n_pcs':None, 'expVar':None})
    
    if expVarThr is not None:
        nC = np.array(np.where(cs>=expVarThr))[0][0]+1
        res['n_pcs'] = nC
        expVar = cs[nC-1]
        res['expVar']= expVar
        
        if figsize is not None:
            plt.hlines(y = expVarThr, xmax=50, xmin=0, colors='green', linestyles='--')
            plt.vlines(x = nC, ymin = 0, ymax= max(cs), colors='green', linestyles='--')
            plt.text(x = nC+0.5, y = 0, s = '{}'.format(nC))
            plt.text(x = 0, y = expVar+0.02, s = '{0:0.3f}'.format(expVar))

    if n_pc is not None:
        expVar = cs[n_pc-1]
        res['expVar']= expVar
        res['n_pcs'] = n_pc
        if figsize is not None:
            plt.hlines(y = expVar, xmax=50, xmin=0, colors='red', linestyles='-.')
            plt.vlines(x = n_pc, ymin = 0, ymax= max(cs), colors='red', linestyles='-.')
            plt.text(x = 0, y = expVar+0.02, s = '{0:0.3f}'.format(expVar))
            plt.text(x = n_pc+0.5, y = 0, s = '{}'.format(n_pc))
        
    if n_pc is not None:
        plt.xlabel('# of PC')
        plt.ylabel('explained variance (cumulative)')
        
    return(res)


def qc_gex_projectlibrary(adata):
    import numpy as np
    import pandas as pd
    import anndata
    
    # Group the data by "project" and "library" using pandas
    grouped = adata.obs.groupby(["libraryID"])
    
    # Create an empty list to store the summary table rows
    summary_rows = []

    # Iterate over each group
    for (libraryID), group_data in grouped:
        
        # Calculate the median values within the group
        chemistry = np.unique(group_data["chemistry"])
        num_cells = group_data["chemistry"].size
        median_counts = np.median(group_data["total_counts"])
        median_features = np.median(group_data["n_genes_by_counts"])
        median_mt_pct = np.median(group_data["pct_counts_mt"])
        median_ribo_pct = np.median(group_data["pct_counts_ribo"])
        median_hb_pct = np.median(group_data["pct_counts_hb"])

        # Create a summary row with the group information and median values
        summary_row = {
            "libraryID": libraryID,
            "chemistry": chemistry,
            "Number of cells": int(num_cells),
            "Median UMI counts": int(median_counts),
            "Median Genes": int(median_features),
            "Median MT %": round(median_mt_pct, 2),
            "Median RB %": round(median_ribo_pct,2),
            "Median HB %": round(median_hb_pct,2)
        }

        # Append the summary row to the list
        summary_rows.append(summary_row)

    # Create a DataFrame from the summary rows
    summary_table = pd.DataFrame(summary_rows)
    return summary_table

def get_HVG(adata, groupby = None, batch_key = None, flavor = 'seurat', min_number_cells = 5, min_disp=0.3, max_disp = 'inf', min_mean=0.01, max_mean=4, n_bins=20):
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scanpy as sc

    if max_disp == 'inf':
        max_disp = np.float16('inf')
        
    HVGdf = pd.DataFrame()
    if groupby is not None:
        listGroups = adata.obs[groupby].unique()

        for g in listGroups:
            print(g)
            adata_g = adata[adata.obs[groupby]==g]
            if adata_g.shape[0] <= min_number_cells:
                print('WARNING: The group {} includes only {} cells. Not considered'.format(g, adata_g.shape[0]))
                continue
            HVGdf_i = sc.pp.highly_variable_genes(adata=adata_g, 
                                    batch_key=batch_key, 
                                    flavor=flavor, 
                                    min_disp=min_disp, 
                                    max_disp=max_disp, 
                                    min_mean=min_mean, 
                                    max_mean=max_mean, 
                                    n_bins=n_bins, inplace=False)

            HVGdf_i = HVGdf_i.add_suffix('_{}'.format(g))
            if batch_key is None:
                HVGdf_i['gene_name'] = adata_g.var_names
                HVGdf_i.set_index('gene_name', inplace=True, drop=True)
                HVGdf_i.index.name = None
            HVGdf = HVGdf.merge(HVGdf_i, how='right', left_index=True, right_index=True)
            

    else:
        HVGdf = sc.pp.highly_variable_genes(adata=adata, 
                                            batch_key=batch_key, 
                                            flavor=flavor, 
                                            min_disp=min_disp, 
                                            max_disp=max_disp, 
                                            min_mean=min_mean, 
                                            max_mean=max_mean, 
                                            n_bins=n_bins,
                                            inplace=False)
        HVGdf['gene_name'] = adata.var_names
        HVGdf.set_index('gene_name', inplace=True, drop=True)
        HVGdf.index.name = None
    
    return(HVGdf)
    

def composition_barplot(adata, xattr, yattr, title, save_pdf, color_dict=None, fig_size=None):
    import pandas as pd
    import matplotlib.pyplot as plt
    
    df = pd.crosstab(adata.obs.loc[:, xattr], adata.obs.loc[:, yattr])
    df = df.div(df.sum(axis=1), axis=0) * 100.0
    if (color_dict is None):
        ax = df.plot(kind = "bar", stacked = True, figsize=fig_size, legend=True, grid=False)
    else:
        ax = df.plot(kind="bar", stacked=True, figsize=fig_size, legend=True, grid=False, color=color_dict)
    ax.figure.subplots_adjust(right=0.9)
    for i in range(len(list(df.index))):
        x_feature = list(df.index)[i]
        feature_count = adata[adata.obs[xattr]==x_feature,:].shape[0]
        ax.annotate(str(feature_count), (ax.patches[i].get_x(), 101))
    ax.annotate("Cell no.", ((ax.patches[i].get_x() +1 ), 101), annotation_clip=False)
    ax.set_title(title,fontsize= 15)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
    plt.xticks(fontsize=12)
    plt.tight_layout()
    plt.show()
    if save_pdf!=None:
        save_pdf.savefig(ax.figure, bbox_inches="tight")
        

def generate_split_dir(cellGroup = None, annotationPath = None, annotation_col_name = None, annotColToCopy = None, targetDir = None, templateNotebook = None):
    import os
    from pyprojroot.here import here

    nb_template = open(templateNotebook, "r") # this line is fixed

    nb_name = os.path.basename(templateNotebook)

    nb_new_path = here('{0}/{1}/{2}'.format(targetDir, cellGroup, nb_name))

    if not os.path.exists(nb_new_path):
        nb_new = open(nb_new_path, "w") # this line change based on the cellGroup label
        line = True
        while line:
            line = nb_template.readline()
            if "cellGroup = 'template'" in line:
                print('{}: cellGroup changed!'.format(nb_name))
                line = line.replace("template", cellGroup)
            if "annotationDFpath = 'annotation_path'" in line:
                print('{}: annotationDFpath changed!'.format(nb_name))
                line = line.replace('annotation_path', annotationPath)
            if "annotationColumns = 'template'" in line:
                print('{}: annotation column changed!'.format(nb_name))
                line = line.replace('template', annotation_col_name)
            if "annColToInclude = None" in line:
                if annotColToCopy is not None:
                    print('{}: annotation to copy changed!'.format(nb_name))
                    line = "\"annColToInclude = '{}'\\n\"".format(annotColToCopy)
            # if "workDir = 'template'" in line:
            #     print('{}: workDir changed!'.format(nb_name))
            #     line = line.replace('template', targetDir+'/'+cellGroup)
            nb_new.write(line)
        nb_template.close()
        nb_new.close()
    else:
        print('WARNING! {} already exist. skipped!'.format(os.path.basename(template_notebook)))



def prepare_data_labelfree_scibmetrics(adata, 
                                embedding_obsm_keys, 
                                pre_integrated_embedding_obsm_key = "X_pca", 
                                neighbor_values = (15, 50, 90), 
                                n_jobs = 10
                               ):
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection 
    from dataclasses import asdict, dataclass
    import os
    import warnings
    from dataclasses import asdict, dataclass
    from enum import Enum
    from functools import partial
    from typing import Any, Callable, Dict, List, Optional, Union
    
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from anndata import AnnData
    from plottable import ColumnDefinition, Table
    from plottable.cmap import normed_cmap
    from plottable.plots import bar
    from sklearn.preprocessing import MinMaxScaler
    from tqdm import tqdm
    
    import scib_metrics
    from scib_metrics.nearest_neighbors import NeighborsOutput, pynndescent   
    
    # (1) Prepare data for comparison
    emb_adatas = {}
    for emb_key in embedding_obsm_keys:
                emb_adatas[emb_key] = AnnData(adata.obsm[emb_key], obs=adata.obs)
                #emb_adatas[emb_key].obs[batch_key] = np.asarray(adata.obs[batch_key].values)
                emb_adatas[emb_key].obsm[pre_integrated_embedding_obsm_key] = adata.obsm[pre_integrated_embedding_obsm_key]
    ## Compute neighbors
    for ad in tqdm(emb_adatas.values(), desc="Computing neighbors"):
        neigh_output = pynndescent(
                ad.X, n_neighbors=max(neighbor_values), random_state=0, n_jobs=n_jobs
            )
        indices, distances = neigh_output.indices, neigh_output.distances
        for n in neighbor_values:
            sp_distances, sp_conns = sc.neighbors._compute_connectivities_umap(
                indices[:, :n], distances[:, :n], ad.n_obs, n_neighbors=n
            )
            ad.obsp[f"{n}_connectivities"] = sp_conns
            ad.obsp[f"{n}_distances"] = sp_distances

    
    return(emb_adatas)

def compare_labelfree_scibmetrics(emb_adatas, 
                                  adata, 
                                batch_key,
                                embedding_obsm_keys, 
                                metric_collection_dict, 
                                pre_integrated_embedding_obsm_key = "X_pca", 
                                neighbor_values = (15, 50, 90), 
                                n_jobs = 10
                               ):
    """
    Function adapted from:
    """
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection 
    from dataclasses import asdict, dataclass
    import os
    import warnings
    from dataclasses import asdict, dataclass
    from enum import Enum
    from functools import partial
    from typing import Any, Callable, Dict, List, Optional, Union
    
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from anndata import AnnData
    from plottable import ColumnDefinition, Table
    from plottable.cmap import normed_cmap
    from plottable.plots import bar
    from sklearn.preprocessing import MinMaxScaler
    from tqdm import tqdm
    
    import scib_metrics
    from scib_metrics.nearest_neighbors import NeighborsOutput, pynndescent

    metric_name_cleaner = {
        "silhouette_label": "Silhouette label",
        "silhouette_batch": "Silhouette batch",
        "isolated_labels": "Isolated labels",
        "nmi_ari_cluster_labels_leiden_nmi": "Leiden NMI",
        "nmi_ari_cluster_labels_leiden_ari": "Leiden ARI",
        "nmi_ari_cluster_labels_kmeans_nmi": "KMeans NMI",
        "nmi_ari_cluster_labels_kmeans_ari": "KMeans ARI",
        "clisi_knn": "cLISI",
        "ilisi_knn": "iLISI",
        "kbet_per_label": "KBET",
        "graph_connectivity": "Graph connectivity",
        "pcr_comparison": "PCR comparison",
    }

    class MetricAnnDataAPI(Enum):
        """Specification of the AnnData API for a metric."""
        #isolated_labels = lambda ad, fn: fn(ad.X, ad.obs[_LABELS], ad.obs[_BATCH])
        #nmi_ari_cluster_labels_leiden = lambda ad, fn: fn(ad.obsp["15_connectivities"], ad.obs[_LABELS])
        #nmi_ari_cluster_labels_kmeans = lambda ad, fn: fn(ad.X, ad.obs[_LABELS])
        #clisi_knn = lambda ad, fn: fn(ad.obsp["90_distances"], ad.obs[_LABELS])
        #silhouette_label = lambda ad, fn: fn(ad.X, ad.obs[_LABELS])
        #graph_connectivity = lambda ad, fn: fn(ad.obsp["15_distances"], ad.obs[_LABELS])
        #silhouette_batch = lambda ad, fn: fn(ad.X, ad.obs[_LABELS], ad.obs[_BATCH])
        pcr_comparison = lambda ad, fn: fn(ad.obsm[pre_integrated_embedding_obsm_key], ad.X, ad.obs[batch_key], categorical=True)
        ilisi_knn = lambda ad, fn: fn(ad.obsp["90_distances"], ad.obs[batch_key])
        #kbet_per_label = lambda ad, fn: fn(ad.obsp["50_connectivities"], ad.obs[_BATCH], ad.obs[_LABELS])
    
    # (1) Add batch_key
    for emb_key in emb_adatas.keys():
                emb_adatas[emb_key].obs[batch_key] = np.asarray(adata.obs[batch_key].values)

    
    # (2) Compare data
    num_metrics = sum(
                [len([key for key, value in asdict(met_col).items() if value is not False]) for met_col in metric_collection_dict.values()]
            )
    results = pd.DataFrame(columns=embedding_obsm_keys + ["Metric Type"])
    for emb_key, ad in tqdm(emb_adatas.items(), desc="Embeddings", position=0, colour="green"):
        pbar = tqdm(total=num_metrics, desc="Metrics", position=1, leave=False, colour="blue")
        for metric_type, metric_collection in metric_collection_dict.items():
            for metric_name, use_metric in asdict(metric_collection).items():
                if use_metric:
                    pbar.set_postfix_str(f"{metric_type}: {metric_name}")
                    metric_fn = getattr(scib_metrics, metric_name)
                    metric_value = getattr(MetricAnnDataAPI, metric_name)(ad, metric_fn)
                    results.loc[metric_name, emb_key] = metric_value
                    results.loc[metric_name, "Metric Type"] = batch_key
                    pbar.update(1)
    return(results)
    
def get_results_labelfree_scibmetrics(df_dicts, 
                                             min_max_scale = True, 
                                             clean_names = True
                               ):

    import pandas as pd
    
    metric_name_cleaner = {
        "silhouette_label": "Silhouette label",
        "silhouette_batch": "Silhouette batch",
        "isolated_labels": "Isolated labels",
        "nmi_ari_cluster_labels_leiden_nmi": "Leiden NMI",
        "nmi_ari_cluster_labels_leiden_ari": "Leiden ARI",
        "nmi_ari_cluster_labels_kmeans_nmi": "KMeans NMI",
        "nmi_ari_cluster_labels_kmeans_ari": "KMeans ARI",
        "clisi_knn": "cLISI",
        "ilisi_knn": "iLISI",
        "kbet_per_label": "KBET",
        "graph_connectivity": "Graph connectivity",
        "pcr_comparison": "PCR comparison",
    }
    
    # Get results
    combined_df = pd.concat(df_dicts.values(), ignore_index=False)
    df = combined_df.transpose()
    df.index.name = "Embedding"
    df = df.loc[df.index != "Metric Type"]
    if min_max_scale:
        # Use sklearn to min max scale
        df = pd.DataFrame(
            MinMaxScaler().fit_transform(df),
            columns=df.columns,
            index=df.index,
        )
    if clean_names:
            df = df.rename(columns=metric_name_cleaner)
    df = df.transpose()
    df["Metric Type"] = combined_df["Metric Type"].values
    per_class_score = df.groupby("Metric Type").mean().transpose()
    per_class_score.columns = per_class_score.columns + ": Agg Score"
    df = df.transpose()
    df.columns = [f'{metric_type}_{col}' for metric_type, col in zip(df.loc['Metric Type'], df.columns)]
    df = pd.concat([df, per_class_score], axis=1)
    df.loc["Metric Type", per_class_score.columns] = [col.replace(": Agg Score", '') for col in per_class_score.columns]
    df = df.transpose().sort_values(by='Metric Type').transpose()
    # Calculate the mean along the columns
    df['Total'] = df.filter(like='Agg Score').drop("Metric Type", axis=0).mean(axis=1)
    df.loc["Metric Type", "Total"] = "Total"
    return(df)

def plot_results_table_labelfree_scibmetrics(df_dicts,
                                             embedding_obsm_keys, 
                                             min_max_scale: bool = True,
                                             show: bool = True, 
                                             save_dir: str = None
                                            ):
    
    import pandas as pd
    import matplotlib.colors as mcolors
    import matplotlib.cm
    from typing import Callable
    from plottable import ColumnDefinition, Table
    from plottable.plots import bar
    import numpy as np
    import matplotlib.pyplot as plt
    
    def custom_normed_cmap(
        s: pd.Series, cmap: mcolors.LinearSegmentedColormap, num_stds: float = 2.5
    ) -> Callable:
        """Returns a normalized colormap function that takes a float as an argument and
        returns an rgba value.
    
        Args:
            s (pd.Series):
                a series of numeric values
            cmap (matplotlib.colors.LinearSegmentedColormap):
                matplotlib Colormap
            num_stds (float, optional):
                vmin and vmax are set to the median ± num_stds.
                Defaults to 2.5.
    
        Returns:
            Callable: Callable that takes a float as an argument and returns an rgba value.
        """
        _median = s.median()
        _std = s.std()
    
        vmin = _median - num_stds * _std
        if vmin < 0:  # Set vmin to 0 if it's negative
            vmin = 0
        vmax = _median + num_stds * _std
    
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    
        return m.to_rgba
    
    num_embeds = len(embedding_obsm_keys)
    cmap_fn = lambda col_data: custom_normed_cmap(col_data, cmap=matplotlib.cm.YlGnBu, num_stds=2.5)
    df = get_results_labelfree_scibmetrics(df_dicts, min_max_scale=min_max_scale)

    # Do not want to plot what kind of metric it is
    plot_df = df.drop("Metric Type", axis=0)
    # Sort by total score
    plot_df = plot_df.sort_values(by="Total", ascending=False).astype(np.float64)
    plot_df["Method"] = plot_df.index
    
    
    # Split columns by metric type, using df as it doesn't have the new method col
    tmp_df = df.copy()
    agg_score_columns = tmp_df.filter(like='Agg Score')
    tmp_df.loc['Metric Type', agg_score_columns.columns] = 'Aggregated score'
    score_cols = tmp_df.columns[tmp_df.loc['Metric Type'].isin(['Total', 'Aggregated score'])]
    other_cols = tmp_df.columns[~tmp_df.loc['Metric Type'].isin(['Total', 'Aggregated score'])]
    column_definitions = [
        ColumnDefinition("Method", width=1.5, textprops={"ha": "left", "weight": "bold"}),
    ]
    
    # Circles for the metric values
    column_definitions += [
        ColumnDefinition(
            col,
            title=col.replace(" ", "\n", 1),
            width=1,
            textprops={
                "ha": "center",
                "bbox": {"boxstyle": "circle", "pad": 0.25},
            },
            cmap=cmap_fn(plot_df[col]),
            group=df.loc["Metric Type", col],
            formatter="{:.3f}",
            border = "left" if i == 0 else None,
        )
        for i, col in enumerate(other_cols)
    ]
    # Bars for the aggregate scores
    column_definitions += [
        ColumnDefinition(
            col,
            width=1,
            title=col.replace(" ", "\n", 1),
            plot_fn=bar,
            plot_kw={
                "cmap": matplotlib.cm.YlGnBu,
                "plot_bg_bar": False,
                "annotate": True,
                "height": 0.9,
                "formatter": "{:.2f}",
            },
            group=df.loc["Metric Type", col],
            border="right" #if i == 0 else None,
        )
        for i, col in enumerate(score_cols)
    ]
    # Allow to manipulate text post-hoc (in illustrator)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, ax = plt.subplots(figsize=(len(df.columns) * 3.3, 5 + 0.3 * num_embeds))
        tab = Table(
            plot_df,
            cell_kw={
                "linewidth": 0,
                "edgecolor": "k",
            },
            column_definitions=column_definitions,
            ax=ax,
            row_dividers=True,
            footer_divider=True,
            textprops={"fontsize": 13, "ha": "center"},
            row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
            col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
            column_border_kw={"linewidth": 1, "linestyle": "-"},
            index_col="Method",
        ).autoset_fontcolors(colnames=plot_df.columns)
    if show:
        plt.show()
    if save_dir is not None:
        fig.savefig(str(save_dir) + "scib_results.svg", facecolor=ax.get_facecolor(), dpi=300)
        
    return tab




################################## functions MLM analysis ##################################

def mean_by_category(acts, category_class):
    import pandas as pd
    
    # Create an empty DataFrame to store the mean values
    mean_by_category = pd.DataFrame(index=acts.var_names)
    # Iterate through each unique category
    for category in acts.obs[category_class].unique():
        # Subset the AnnData object for the current category
        subset_acts = acts[acts.obs[category_class] == category]
        # Calculate the mean expression for each gene within the category
        mean_expression = subset_acts.X.mean(axis=0)
        # Add the mean expression values to the DataFrame
        mean_by_category[category] = list(mean_expression)
    return mean_by_category

def filter_low_represented_cell_group(adata=None, min_nCell = 20, annotation_col = 'Level2', disease_col = 'disease', reference_group = 'healthy'):
    r = adata.obs.groupby([annotation_col, disease_col]).size().reset_index()
    r.columns = [annotation_col,disease_col,'size']
    rfilt = r.query("size >= @min_nCell")
    cg_reference = rfilt[annotation_col][rfilt[disease_col] == reference_group].unique()
    finalFilt = rfilt[rfilt[annotation_col].isin(cg_reference)]

    print("filtering")
    adataFilt = adata[adata.obs[annotation_col].isin(finalFilt[annotation_col])].copy()
    
    print("The following cell groups have been discarded")
    [print(c) for c in adata.obs[annotation_col][~adata.obs[annotation_col].isin(finalFilt[annotation_col])].unique()]
    print(f"We kept {finalFilt[annotation_col].unique().shape[0]} cell groups, removing a total of {adata.shape[0] - adataFilt.shape[0]} cells")
    
    return adataFilt
    

def RelativeDiff_mean_by_category(mean_by_category, reference):
    import pandas as pd
    import numpy as np

    # Create an empty DataFrame to store the results
    mean_by_category_RD = pd.DataFrame(index=mean_by_category.index, columns=mean_by_category.columns)
    mean_by_category_RD_subset = pd.DataFrame(index=mean_by_category.index, columns=mean_by_category.columns)

    # Iterate over rows and columns to calculate the maximum value for each comparison
    for row in mean_by_category.index: # pathway
        for col in mean_by_category.columns: # diseases
            if col != reference:
                mean_by_category_RD.loc[row, col] = (( mean_by_category.loc[row, col] - mean_by_category.loc[row, reference] ) / abs(mean_by_category.loc[row, reference]))*100
                if (mean_by_category.loc[row, col] < 0) & (mean_by_category.loc[row, reference] < 0):
                    mean_by_category_RD_subset.loc[row, col] = np.nan
                else:
                    mean_by_category_RD_subset.loc[row, col] = -1

    # Drop the 'healthy' column from the result DataFrame
    mean_by_category_RD = mean_by_category_RD.drop(reference, axis=1)
    mean_by_category_RD_subset = mean_by_category_RD_subset.drop(reference, axis=1)
    # You can choose another data type if needed
    mean_by_category_RD = mean_by_category_RD.astype(np.float64)  
    mean_by_category_RD_subset = mean_by_category_RD_subset.astype(np.float64)

    # Apply the log10 transformation to all values in the DataFrame
    mean_by_category_RD_log10 = (mean_by_category_RD).applymap(lambda x: -1 * np.log10(1 + np.abs(x)) if x < 0 else np.log10(1 + np.abs(x)))

    return mean_by_category_RD_log10, mean_by_category_RD, mean_by_category_RD_subset


def mscatter(x,y,ax=None, m=None, **kw):    
    import matplotlib.markers as mmarkers
    
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc
    
################################################# PATIENT CLASSIFIER #####################################################

def aggregating_features(Z = None, obsDF = None, mode = 'mean', obs_names_col = [], min_observation = 0):

    import pandas as pd
    import scanpy as sc
    import numpy as np
    
    Zdf = pd.DataFrame(Z)
    for c in obsDF.columns:
        Zdf[c] = obsDF[c].tolist()

    grpDF = Zdf.groupby(obsDF.columns.tolist(), observed = True)

    nCount = grpDF.size().to_frame('n_observation')
    
    if mode in ['mean','avarage']:
        Zaggr = grpDF.mean()
    elif mode == 'sum':
        Zaggr = grpDF.sum()
    else:
        raise ValueError(f"mode {mode} not supported. Available mode are 'mean' or 'sum'")

    grpObs = pd.DataFrame(Zaggr.index.tolist(), columns=obsDF.columns.tolist()).merge(pd.DataFrame(nCount).reset_index(), on = obsDF.columns.tolist())

    if len(obs_names_col) == 0:
        grpAdata  = sc.AnnData(X = np.array(Zaggr), obs = grpObs)
    elif all([c in obsDF.columns.tolist() for c in obs_names_col]):
        grpObs.index = grpObs[obs_names_col].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        grpAdata  = sc.AnnData(X = np.array(Zaggr), obs = grpObs)
    else:
        raise ValueError(f"Impossible to use {obs_names_col} as index. It's not present in obsDF")

    if min_observation > 0:
        grpAdata = grpAdata[grpAdata.obs.n_observation >= min_observation]
    return grpAdata

def train_patient_classifier(adataTrain = None, cell_type_col = None, y_true_col = None, max_iter=10000,random_state = 25, model = 'LinearSVC', kargs_model = None):
    
    import pandas as pd
    from tqdm import tqdm
    from sklearn.svm import LinearSVC, SVC
    from sklearn.neighbors import KNeighborsClassifier

    from sklearn.metrics import balanced_accuracy_score
    
    clfList = dict()
    trainAccuracy = []
    for ct_ in tqdm(adataTrain.obs[cell_type_col].unique()):
        clfList[ct_] = dict()
        X_i = adataTrain.X[adataTrain.obs[cell_type_col] == ct_]
        y_true_i = adataTrain.obs[y_true_col][adataTrain.obs[cell_type_col] == ct_]

        if model == 'LinearSVC':
            clfList[ct_]['clf'] = LinearSVC(**kargs_model).fit(X_i, y_true_i) # max_iter=max_iter, dual = True, random_state = random_state
        elif model == 'SVC':
            clfList[ct_]['clf'] = SVC(**kargs_model).fit(X_i, y_true_i) # max_iter=max_iter, random_state = random_state
        elif model == 'KNeighborsClassifier':
            clfList[ct_]['clf'] = KNeighborsClassifier(**kargs_model).fit(X_i, y_true_i) # n_neighbors = 5, weights='distance', n_jobs = -1
            
            
        clfList[ct_]['bAcc'] = balanced_accuracy_score(y_true = y_true_i, y_pred = clfList[ct_]['clf'].predict(X_i))
        clfList[ct_]['nObs'] = len(y_true_i)

    return clfList
    
def vote_patient_disease(adataTest = None, clfList = None, cell_type_col = None, sample_id_col = None):

    from tqdm import tqdm
    import numpy as np
    import pandas as pd
    
    classificationDF = pd.DataFrame()
    for ct_ in tqdm(adataTest.obs[cell_type_col].unique()):
        X_i = adataTest.X[adataTest.obs[cell_type_col] == ct_]
        PID_i = adataTest.obs[sample_id_col][adataTest.obs[cell_type_col] == ct_]
        if ct_ not in clfList:
            print(f"{ct_} is missing in training set")
            continue
        DF_i = pd.DataFrame.from_dict({
            sample_id_col: PID_i,
            f"{ct_}_prediction": clfList[ct_]['clf'].predict(X_i)
        })
        if classificationDF.shape[0] == 0:
            classificationDF = DF_i
        else:
            classificationDF = classificationDF.merge(DF_i, how='outer', on = sample_id_col)
    classificationDF['firstChoice'] = ''
    classificationDF['firstChoice_perc'] = np.nan
    classificationDF['secondChoice'] = ''
    classificationDF['secondChoice_perc'] = np.nan
    
    for i in tqdm(range(classificationDF.shape[0])):
        vote_i = classificationDF.loc[i,classificationDF.columns !=sample_id_col].value_counts()
        vote_i /= (vote_i.sum() / 100)
        res_i = vote_i.sort_values(ascending=False)
        classificationDF.loc[i,'firstChoice'] = res_i.index[0]
        classificationDF.loc[i,'firstChoice_perc'] = res_i.iloc[0]
        if res_i.shape[0] > 1:
            classificationDF.loc[i,'secondChoice'] = res_i.index[1]
            classificationDF.loc[i,'secondChoice_perc'] = res_i.iloc[1]

    return classificationDF