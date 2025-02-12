import pandas as pd
import scanpy as sc
import os
import sys
import numpy as np
import sys
import tqdm
import math

import SEACells

# ---------------------------------------------------------------------------------------------------
# CUSTOM FUNCTIONS
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

import copy
def get_soft_assignments_modified(self):
    """Compute soft SEACells assignment."""
    archetype_labels = self.get_hard_archetypes()
    A = copy.deepcopy(self.A_.T)

    labels = []
    weights = []
    for _i in range(5):
        l = A.argmax(1)
        #labels.append(archetype_labels[l])
        seacell_name = [f"SEACell-{i}" for i in l]
        labels.append(seacell_name)
        weights.append(A[np.arange(A.shape[0]), l])
        A[np.arange(A.shape[0]), l] = -1

    weights = np.vstack(weights).T
    labels = np.vstack(labels).T

    soft_labels = pd.DataFrame(labels)
    soft_labels.index = self.ad.obs_names

    return soft_labels, weights

# ---------------------------------------------------------------------------------------------------
#                            IMPORT PARAMETERS
# ---------------------------------------------------------------------------------------------------
# getting input parameters
sampleAdataFN = sys.argv[1]
SCxSEACell = int(sys.argv[2])

# ---------------------------------------------------------------------------------------------------
#                            LOAD DATA
# ---------------------------------------------------------------------------------------------------
# --- LOAD DATA --- 
np.random.seed(42)
# reading sample adata
sampleAdata = sc.read_h5ad(sampleAdataFN)

sample_id = sampleAdata.obs.sampleID.unique()
assert(sample_id.shape[0] == 1)
sample_id = sample_id[0]
print(f"******* Computing SEACell for {sample_id} *******")

# ---------------------------------------------------------------------------------------------------
#                            PROCESS ADATA
# ---------------------------------------------------------------------------------------------------


# --- PREPARE OBJECT FOR SEACell --- 
sampleAdata.X = sampleAdata.X.astype(float)
sc.pp.normalize_per_cell(sampleAdata)
sc.pp.log1p(sampleAdata)

# ---- SELECT HVG ----
## Parameters 
HVG_parameters = {
    "flavor":'seurat', 
    "n_top_genes":3000
}

## Compute HVG
if sampleAdata.obs.libraryID.unique().shape[0] > 1: 
    batch_key = "libraryID"
else: 
    batch_key = None
sc.pp.highly_variable_genes(adata=sampleAdata, 
                                    batch_key=batch_key, 
                                    flavor=HVG_parameters["flavor"], 
                                    n_top_genes = HVG_parameters["n_top_genes"])

# ---- PCA ----
sc.tl.pca(sampleAdata, svd_solver="arpack", n_comps=50, use_highly_variable=True)
minExpVar = 0.9
n_pcs = cumulative_explained_variance(sampleAdata, expVarThr = minExpVar, n_pc=None, figsize=None)['n_pcs']
X_pca_npcs = sampleAdata.obsm['X_pca'][:, :n_pcs]
sampleAdata.obsm['X_pca_npcs'] = X_pca_npcs

# ---------------------------------------------------------------------------------------------------
#                            RUN SEACell
# ---------------------------------------------------------------------------------------------------

# ---- SEACell parameters ----
n_SEACells = math.ceil(sampleAdata.n_obs / SCxSEACell)

print(f"Will be generated {n_SEACells} SEAcells")

build_kernel_on = 'X_pca_npcs'                          
n_waypoint_eigs = 7 

# ---- SEACell ----
model = SEACells.core.SEACells(sampleAdata, 
          n_neighbors = 10, #default is 15, changed to solve an error
          build_kernel_on=build_kernel_on, 
          n_SEACells=n_SEACells, 
          n_waypoint_eigs=n_waypoint_eigs,
          convergence_epsilon = 1e-5)
model.construct_kernel_matrix()
M = model.kernel_matrix
model.initialize_archetypes()
model.fit(min_iter=10, max_iter=10000)
# ---- Import results ----
labels,weights = get_soft_assignments_modified(model)
max_values = np.max(weights, axis=1).tolist()
non_trivial_assig = (model.A_.T > 0.1).sum(axis=1)
sampleAdata.obs['SEACell'] = sample_id + "_" + sampleAdata.obs['SEACell']
sampleAdata.obs['SEACell_weight'] = max_values
sampleAdata.obs['SEACell_non_trivial_assig'] = non_trivial_assig
# ---- Generate DF ----
SEACell_DF = sampleAdata.obs[['SEACell', 'SEACell_weight', 'SEACell_non_trivial_assig']]
SEACell_DF.index = sampleAdata.obs.index
# ---------------------------------------------------------------------------------------------------
#                            ADD QC results
# ---------------------------------------------------------------------------------------------------
compactness = SEACells.evaluate.compactness(sampleAdata, build_kernel_on)
separation = SEACells.evaluate.separation(sampleAdata, build_kernel_on, nth_nbr=1)
SEACell_DF = SEACell_DF.merge(compactness, on='SEACell', how='left')
SEACell_DF = SEACell_DF.merge(separation, on='SEACell', how='left')
SEACell_DF.index = sampleAdata.obs.index
# ---------------------------------------------------------------------------------------------------
#                            SAVE SEACell results
# ---------------------------------------------------------------------------------------------------
SEACell_DF.to_pickle(f"{sample_id}.seacell.pkl")