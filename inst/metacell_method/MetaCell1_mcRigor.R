

library(Seurat)
#library(doParallel)
library(metacell)




MetaCell1_mcRigor <- function(X,  # gene by cell count matrix
                              cell.metadata = NULL,
                              mc_dir = 'MC1/',
                              mc_id = 'mc1',
                              GENE_USE = NULL,
                              EXCLUDED_GENE_NAMES = NULL,
                              EXCLUDED_GENE_PATTERNS = NULL,
                              CELL_MIN_UMI = 1,
                              k_knn = 50,
                              initial_knn_amp = 6,
                              outlier_detect = T,
                              seed = 123456,
                              filter_gene = T,
                              n_downsamp_gstat = NULL,
                              do_downsamp = T,
                              ...){
  
  N.c <- ncol(X)
  
  if(is.null(colnames(X))){
    warning("colnames(X) is Null, \nGene expression matrix X is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' ")
    colnames(X) <- paste("cell", 1:N.c, sep = "_")
  }
  
  cell.ids <- colnames(X)
  
  
  ###
  
  if(!dir.exists(mc_dir)) dir.create(mc_dir)
  scdb_init(mc_dir, force_reinit=T)
  
  scdb_add_mat(id = mc_id, mat = scm_new_matrix(X, cell_metadata = cell.metadata))
  
  if(!dir.exists(paste0(mc_dir, "figs/"))) dir.create(paste0(mc_dir, "figs/"))
  scfigs_init(paste0(mc_dir, "figs/"))
  
  # ignore some bad genes
  mat = scdb_mat(mc_id)
  nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
  
  if (!is.null(GENE_USE)) {
    bad_genes = setdiff(nms, GENE_USE)
  } else {
  ig_genes = c(grep("^IGJ", nms, v=T), 
               grep("^IGH",nms,v=T),
               grep("^IGK", nms, v=T), 
               grep("^IGL", nms, v=T))
  
  bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),
                       "NEAT1","TMSB4X", "TMSB10", ig_genes))
  
  if (!is.null(EXCLUDED_GENE_NAMES)) bad_genes = unique(c(bad_genes, EXCLUDED_GENE_NAMES))
  if (!is.null(EXCLUDED_GENE_PATTERNS)) bad_genes = unique(c(bad_genes, 
                                                             do.call(c, lapply(EXCLUDED_GENE_PATTERNS, function(gp) grep(gp, nms, v=T)))))
  }
  
  mcell_mat_ignore_genes(new_mat_id=mc_id, mat_id=mc_id, bad_genes, reverse=F) 
  
  mcell_mat_ignore_small_cells(mc_id, mc_id, CELL_MIN_UMI)
  
  if (filter_gene) {
    
    if (!is.null(n_downsamp_gstat)) tgconfig::set_param('scm_n_downsamp_gstat', n_downsamp_gstat, package = 'metacell')
    mcell_add_gene_stat(gstat_id=mc_id, mat_id=mc_id, force=T)
    if (!is.null(n_downsamp_gstat)) tgconfig::set_param('scm_n_downsamp_gstat', NULL, package = 'metacell')
    
    mcell_gset_filter_varmean(gset_id=mc_id, gstat_id=mc_id, T_vm=0.08, force_new=T)
    mcell_gset_filter_cov(gset_id = mc_id, gstat_id=mc_id, T_tot=100, T_top3=2)
    
  } else {
    
    mat = scdb_mat(mc_id)
    gset = rep(1, length(mat@genes))
    names(gset) = mat@genes
    scdb_add_gset(mc_id, gset = gset_new_gset(gset, 'all'))
    
  }
  
  mcell_add_cgraph_from_mat_bknn(mat_id=mc_id, 
                                 gset_id =mc_id, 
                                 graph_id=mc_id,
                                 K=k_knn * initial_knn_amp,
                                 dsamp=do_downsamp)
  
  mcell_coclust_from_graph_resamp(
    coc_id=mc_id, 
    graph_id=mc_id,
    min_mc_size=1, 
    p_resamp=0.75, n_resamp=500)
  
  mcell_mc_from_coclust_balanced(
    coc_id=mc_id, 
    mat_id= mc_id,
    mc_id= mc_id, 
    K=k_knn, min_mc_size=1, alpha=2)
  
  if (outlier_detect) {
    mcell_mc_split_filt(new_mc_id=mc_id, 
                        mc_id=mc_id, 
                        mat_id=mc_id,
                        T_lfc=3, plot_mats=F)
  }
  
  mc = scdb_mc(mc_id)
  
  sc_membership = rep(NA, N.c)
  names(sc_membership) = cell.ids
  
  sc_membership[names(mc@mc)] = mc@mc

  gamma = round(length(mc@mc) / length(table(sc_membership[!is.na(sc_membership)])))
  
  res = list(gamma = gamma, k_knn = k_knn, sc_membership = sc_membership)
  
  return(res)
}




