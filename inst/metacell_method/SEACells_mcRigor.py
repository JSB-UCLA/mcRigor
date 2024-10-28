
#%%
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import scanpy as sc
import scvelo as scv


#%%
def SEACells_mcRigor(input_file,
                     output_file = 'seacells_cell_membership_rna.csv',
             prePro = False,
             Gamma = range(100, 10, -1),
             reduction_key = "X_pca",
             dim_str = "1:50",
             n_features = 2000,
             annotations = None,
             min_metacells = 1,
             n_waypoint_eigs = 10, # Number of eigenvalues to consider when initializing metacells
             min_iter=10,
             max_iter=100,
             k_knn = 30
             ):
    

    print('input is "', input_file)
    print('Gamma include "', [i for i in Gamma])
    print('dims are "', dim_str)
    print('reduction_key is"', reduction_key)
    
    if input_file.endswith(".h5ad"):
        adata = sc.read_h5ad(input_file)
        
        if adata.raw is not None:
            adata.X = adata.raw  # we only load raw counts, We always normalize .X prior to compute PCA if prePro is asked or reduction_key absent  
            del adata.raw

    else:
        adata = sc.read_csv(input_file + 'counts.csv')
        metadata = pd.read_csv(input_file + 'metadata.csv', index_col = 0)
        adata.obs = metadata

        from scipy import sparse
        sparse_X = sparse.csr_matrix(adata.X)
        adata.X = sparse_X

        if (reduction_key == 'X_svd' and reduction_key not in adata.obsm.keys()):
            reduction = pd.read_csv(input_file + 'reduction.csv', index_col = 0)
            adata.obsm['X_svd'] = reduction.values

    # The dtype of X is no longer set to float32 in scampy. 
    # While anndata2ri produces float64, the majority of h5ad objects available online are float32.
    # We choose to set the type to float32
    
    adata.X = adata.X.astype("float32")

    ####################
    
    dim_str_list = dim_str.split(":") # range input is interesting when using SEACells for ATAC data for which 1 component is often discarded
    # 1 to given components are used when only one number is given to dim
    if (len(dim_str_list)<2):
        dim_str_list += "1"
        dim_str_list.reverse()
    
    # Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
    # This step should be performed after filtering
    raw_ad = sc.AnnData(adata.X)
    raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
    adata.raw = raw_ad

    if annotations is None:
        annotations = "SEACell_batch"
        adata.obs["SEACell_batch"] = "allcells"

    cell_membership = pd.DataFrame(index=adata.obs.index)
    for gamma in Gamma:
        try:
            print("Identify metacells using SEACells. gamma = " + str(gamma))  

            for anno in adata.obs[annotations].unique():
                adata_label = adata[adata.obs[annotations] == anno,]
                n_SEACells = round(len(adata_label)/gamma)
                
                if n_SEACells < min_metacells:
                    n_SEACells = min_metacells
                    
                if n_SEACells == 1:
                    adata_label.obs['SEACell'] = anno + "-" + "SEACell-1"
                    if anno == adata.obs[annotations].unique()[0]:
                        seacells_res = adata_label.obs["SEACell"]
                    else:
                        seacells_res = pd.concat([seacells_res,adata_label.obs["SEACell"]])
                    continue
                
                print("Identify "+ str(n_SEACells) + " metacells using SEACells...")

                if (reduction_key == 'X_svd'):
                    prePro = False
                    if (reduction_key not in adata_label.obsm.keys()):
                        raise Exception("No SVD reduction for scATAC data!")
                elif (prePro or reduction_key not in adata_label.obsm.keys()):
                    print("Preprocess the data...")
                    print("Normalize cells and compute highly variable genes...")
                    sc.pp.normalize_per_cell(adata_label)
                    sc.pp.log1p(adata_label)
                    sc.pp.highly_variable_genes(adata_label, n_top_genes=n_features)
                
                    print("Compute principal components")
                    sc.tl.pca(adata_label, n_comps=int(dim_str_list[1]), use_highly_variable=True)
                    reduction_key = "X_pca"
                
                build_kernel_on = reduction_key
                dim_lwb = int(dim_str_list[0])-1  # The first dimension of SVD for scATAC should be removed beforehand
                dim_upb = min(int(dim_str_list[1]), adata_label.obsm[build_kernel_on].shape[1])
                adata_label.obsm[build_kernel_on] = adata_label.obsm[build_kernel_on][:,range(dim_lwb, dim_upb)]
                
            
                min_metacells = min(min_metacells,adata_label.n_obs)
            
            
                if n_SEACells < n_waypoint_eigs:
                    n_waypoint_eigs_label = n_SEACells
                else:
                    n_waypoint_eigs_label = n_waypoint_eigs
            
            
                model = SEACells.core.SEACells(adata_label, 
                build_kernel_on=build_kernel_on, 
                n_SEACells=n_SEACells, 
                n_neighbors = k_knn,
                n_waypoint_eigs=n_waypoint_eigs_label,
                convergence_epsilon = 1e-5)
                
                model.construct_kernel_matrix()
                # M = model.kernel_matrix
                model.initialize_archetypes()
                #SEACells.plot.plot_initialization(ad, model)
                model.fit(min_iter=min_iter, max_iter=max_iter)
                #model.plot_convergence()
                adata_label.obs['SEACell'] = "mc" + str(gamma) +"-"+ anno + "-" + adata_label.obs['SEACell']  
                
                if anno == adata.obs[annotations].unique()[0]:
                    seacells_res = adata_label.obs["SEACell"]
                else:
                    seacells_res = pd.concat([seacells_res,adata_label.obs["SEACell"]])
                

            adata.obs['SEACell'] = seacells_res.reindex(adata.obs_names)

            cell_membership[str(gamma)] = adata.obs['SEACell']
            cell_membership.to_csv(output_file)

        except:
            print('Something went wrong for gamma = ' + str(gamma))
        
    print("All finished!")




# # %%
# input_file = '/home/pan/wrapup/cd34_multiome_rna.h5ad'
# SEACells_mcRigor(input_file = input_file, Gamma = range(13,10,-1), 
#                  output_file = 'seacells_cell_membership_rna_2.csv',
#                  reduction_key = "X_pca")

# # %%
# input_file = '/home/pan/wrapup/cd34/input/cd34_multiome_atac.h5ad'
# # adata = sc.read_h5ad(input_file)
# SEACells_mcRigor(input_file = input_file, Gamma = range(200,100,-1), 
#                  output_file = 'seacells_cell_membership_atac_1.csv',
#                  reduction_key='X_svd')

# %%
input_file = '/home/pan/wrapup/covid_pbmc/input/covid/'
SEACells_mcRigor(input_file = input_file, Gamma = range(100,1,-1), 
                 output_file = '/home/pan/wrapup/covid_pbmc/seacells_cell_membership_rna_covid.csv',
                 reduction_key='X_pca')

# %%
input_file = '/home/pan/wrapup/covid_pbmc/input/healthy/'
SEACells_mcRigor(input_file = input_file, Gamma = range(100,1,-1), 
                 output_file = '/home/pan/wrapup/covid_pbmc/seacells_cell_membership_rna_healthy.csv',
                 reduction_key='X_pca')


