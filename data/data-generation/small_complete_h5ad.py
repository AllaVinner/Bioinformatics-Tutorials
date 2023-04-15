import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


def generate(num_cells = 100, num_genes = 20000):

    counts = csr_matrix(np.random.poisson(1, size=(num_genes, num_cells)), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

    cell_types = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
    adata.obs["cell_type"] = pd.Categorical(cell_types)







