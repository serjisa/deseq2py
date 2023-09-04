import scanpy as sc
import numpy as np
import pandas as pd

from .. import logging as logg

def pseudobulk(adata, batch_key, keep_columns=None, layer=None, use_raw=False, sparse=False):
    if use_raw:
        if not (layer is None):
            logg.warn("Can't combine both `layer` and `use_raw`, using only `adata.raw.X`")
        matrix = adata.raw.X
        var_names = list(adata.raw.var_names)
    elif layer is None:
        matrix = adata.X
        var_names = list(adata.var_names)
    else:
        matrix = adata.layers[layer]
        var_names = list(adata.var_names)

    samples = list(set(adata.obs[batch_key]))
    adata_pseudobulk = sc.AnnData(
        X=np.array([
            matrix[list(adata.obs[batch_key] == sample)].sum(axis=0).A[0]
            for sample in samples
        ]),
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame(index=samples),
    )
    
    if not (keep_columns is None):
        for column in keep_columns:
            if len(np.unique(adata.obs[[batch_key, column]], axis=0)) != len(samples):
                logg.warn(f"More than one label from adata.obs[{column}] corresonds to some single batch, skipping...")
            else:
                batch_to_column = dict(adata.obs[[batch_key, column]].astype(str).values)
                adata_pseudobulk.obs[column] = [batch_to_column[sample] for sample in adata_pseudobulk.obs_names]
              
    if sparse:
        from scipy.sparse import csr_matrix
      
        adata_pseudobulk.X = csr_matrix(adata_pseudobulk.X)
        
    return adata_pseudobulk
