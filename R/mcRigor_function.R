




#' @title mcRigor_DETECT
#'
#' @description Main functionality 1: To detect dubious metacells for a given metacell partition
#'
#'
#' @param obj_singlecell Seurat object of the original single cells.
#' @param cell_membership A dataframe, each column of which contains the metacell membership of single cells under a given gamma (granularity level). 
#' The column names should be the corresponding gamma's (in the character type). The row names should be the single cell id's.
#' @param tgamma The target gamma value --- the gamma (in the character type) corresponding to the metacell partition where dubious metacells detection will be performed.
#' tgamma has to be in the column names of cell_membership. If tgamma = NULL, dubious metacells detection will be performed on the metacell partition represented by the first column of cell_membership.
#' @param assay_type The type of data assay yuo are using, depending on which different normalization would be used.
#' @param aggregate_method The method to aggregate single cell profiles into metacell profiles
#' @param output_file The directory for output saving
#' @param Nrep The number of permutation repetitions we use when deriving the null.
#' @param test_cutoff The test size for dubious metacell detection testing
#' @param thre_smooth A boolean indicating whether to smooth the threshold function
#' @param thre_bw If thre_smooth is True, thre_bw specifies the bandwidth for smoothing
#' @param cor_method The method for gene correlation calculation description
#' @param gene_filter A proportion. Genes expressed lower than this proportion will be filtered out.
#' @param feature_use The number of genes to use in metacell testing.
#' @param prePro A boolean indicating whether to normalize obj_singlecell for preprocessing.
#' @param check_purity A boolean indicating whether to calculate the metacell purity of specific fields or not. 
#' @param fields A vector of the fields of interest, ex. celltype. It should be a subset of obj_singlecell's meta.data.
#' @param draw  A boolean indicating whether to visualize the mcRigor results
#' @param pur_metric Can be NULL or a metadata variable name, ex. cell type.
#' @param step_save A boolean indicating whether to save the outputs step by step
#' 
#' 
#' @return A list containing the following fields: 
#' \item{best_graining_level}{The optimal graining level selected by mcRigor}
#' \item{opt_metacell}{The metacell object build under the optimal graining level}
#' \item{thre}{The thresholds for dubious metacell detection}
#' \item{TabMC}{A dataframe containing the permutation results, elements to calculate the test statistics mcDiv and mcDiv null}
#' \item{test_plot}{The scatter plots demonstrating the mcDiv values and the obtained thresholds for dubious metacell detection}
#' \item{purity_plot}{A violin plot showing the purity distribution of the pur_metric covariate in dubious metacells and trustworthy metacells}
#'


#
mcRigor_DETECT <- function(obj_singlecell, 
                           cell_membership = NULL, 
                           tgamma = NULL,
                           assay_type = c('RNA', 'ATAC'),
                           aggregate_method = c('mean', 'sum', 'geom'),
                           output_file = NULL,
                           Nrep = 1,
                           gene_filter = 0.1, 
                           feature_use = 2000,
                           cor_method = c('pearson', 'spearman'),
                           prePro = T,
                           test_cutoff = 0.01, 
                           thre_smooth = T, 
                           thre_bw = 1/6,
                           draw = T, 
                           pur_metric = NULL,
                           check_purity = T,
                           fields = NULL,
                           step_save = T){
  
  assay_type = match.arg(assay_type)
  cor_method = match.arg(cor_method)
  aggregate_method = match.arg(aggregate_method)
  
  if (prePro){
    cat("Normalizing data...\n")
    if (assay_type == 'ATAC') obj_singlecell = Signac::RunTFIDF(obj_singlecell, verbose = F)
    else obj_singlecell = Seurat::NormalizeData(obj_singlecell, verbose = F)
  }
  
  if (is.null(fields)) {
    fields = sapply(colnames(obj_singlecell@meta.data),
                    function(X) {is.character(obj_singlecell[[X]][,1]) | is.factor(obj_singlecell[[X]][,1])})
    fields = names(fields[fields])
  }
  
  if (assay_type == 'ATAC'){
    feature_cutoff = paste0('q', floor((1 - feature_use / dim(obj_singlecell)[1]) * 100))
    obj_singlecell = Signac::FindTopFeatures(obj_singlecell, min.cutoff = feature_cutoff, verbose = F)
    top_gene = Seurat::VariableFeatures(obj_singlecell)
  } else {
    obj_singlecell = Seurat::FindVariableFeatures(obj_singlecell, nfeatures = feature_use, verbose = F)
    top_gene = Seurat::VariableFeatures(obj_singlecell)
  }
  
  if (is.null(tgamma)) tgamma = colnames(cell_membership)[1]
  
  Gammas = colnames(cell_membership)
  
  TabMC_list = vector(mode = 'list', length = length(Gammas))
  names(TabMC_list) = Gammas
  
  for (gamma in Gammas) {
    
    named_membership = cell_membership[[gamma]]
    names(named_membership) = rownames(cell_membership)
    assigned_cellid = names(named_membership)[named_membership != '' & !is.na(named_membership)]
    named_membership = named_membership[assigned_cellid]
    
    obj_metacell = mcRigor_buildmc(obj_singlecell[, assigned_cellid], sc_membership = named_membership,
                                   assay_type = assay_type, fields = fields, aggregate_method = 'sum',
                                   doNorm = F)
    obj_metacell[['ZeroRate']] = 1 - obj_metacell[[paste0('nFeature_', Seurat::DefaultAssay(obj_metacell))]] / nrow(obj_metacell)
    
    ##########
    
    metacell_id = colnames(obj_metacell)
    
    mc_stats = data.frame(matrix(NA, nrow = length(metacell_id), ncol = 3+2*Nrep))
    rownames(mc_stats) = metacell_id
    colnames(mc_stats) = c('size', 'T_org', 'T_colperm',
                           paste0('T_rowperm', 1:Nrep), paste0('T_bothperm', 1:Nrep))
    #
    cat('gamma =', gamma, '\n')
    pb <- utils::txtProgressBar(style=3)
    pbid = 0
    
    for (select.label in metacell_id) {
      
      select_sc = names(named_membership)[named_membership == select.label]
      mc_stats[select.label, 'size'] = length(select_sc)
      
      if (length(select_sc) == 1) next
      
      if(dim(obj_singlecell)[1] > dim(obj_singlecell)[2]){
        counts = Seurat::GetAssayData(object = obj_singlecell, layer = 'data')[top_gene, select_sc]
      } else {
        select_mc = obj_singlecell[, select_sc]
        counts = Seurat::GetAssayData(select_mc, layer = 'data')
      }
      
      counts = as.matrix(counts)
      # counts = apply(counts, 2, function(x) x/sqrt(sum(x^2)))
      
      counts_var = counts[apply(counts, 1, 
                                function(x) length(which(x>0)) > dim(counts)[2] * gene_filter), ]
      
      if (cor_method == 'pearson') dat = scale_cpp(t(counts_var))
      else {
        dat = apply(t(counts_var), 2, rank)
        dat = scale_cpp(dat)
      }
      
      p = dim(dat)[2]
      n = dim(dat)[1]
      
      mc_stats[select.label, 'T_org'] = mc_indpd_stats_cpp(dat)
      mc_stats[select.label, 'T_colperm'] = mc_indpd_stats_cpp(colwise_perm_cpp(dat))
      
      T_perm = sapply(1:Nrep, function(i) {
        x = rowwise_perm_cpp(dat)
        c(mc_indpd_stats_cpp(x), mc_indpd_stats_cpp(colwise_perm_cpp(x)))
      })
      mc_stats[select.label, paste0('T_rowperm', 1:Nrep)] = T_perm[1,]
      mc_stats[select.label, paste0('T_bothperm', 1:Nrep)] = T_perm[2,]
      
      pbid = pbid + 1
      utils::setTxtProgressBar(pb, pbid / length(metacell_id))
      
    }
    
    close(pb)
    
    purity_fields = colnames(obj_metacell@meta.data)[grep('_purity', colnames(obj_metacell@meta.data))]
    TabMC_list[[gamma]] =  cbind(as.numeric(gamma), 
                                 obj_metacell[['ZeroRate']],
                                 obj_metacell[[purity_fields]], 
                                 mc_stats)
    rownames(TabMC_list[[gamma]]) = rownames(mc_stats)
    
    if (step_save){
      TabMC = do.call(rbind, unname(TabMC_list))
      colnames(TabMC)[1] = 'gamma'
      colnames(TabMC)[2] = 'ZeroRate'
      
      TabMC$TT_div = TabMC$T_org / TabMC$T_colperm
      TabMC$TT_div[is.na(TabMC$TT_div)] = 1
      
      if (is.null(output_file)) 
        output_file = paste0('Tabmc_', assay_type, '_detect.RData')
      save(TabMC, file = output_file)
    }
    
  }
  
  TabMC = do.call(rbind, unname(TabMC_list))
  colnames(TabMC)[1] = 'gamma'
  colnames(TabMC)[2] = 'ZeroRate'
  
  TabMC$TT_div = TabMC$T_org / TabMC$T_colperm
  TabMC$TT_div[is.na(TabMC$TT_div)] = 1
  
  if (is.null(output_file)) 
    output_file = paste0('Tabmc_', assay_type, '_detect.RData')
  save(TabMC, file = output_file)
  
  ###########
  
  test_res = mcRigor_threshold(TabMC,
                               test_cutoff = test_cutoff, 
                               thre_smooth = thre_smooth, thre_bw = thre_bw, 
                               draw = draw, pur_metric = pur_metric)
  Thre = test_res$threshold
  
  ##
  
  named_membership = cell_membership[[as.character(tgamma)]]
  names(named_membership) = rownames(cell_membership)
  assigned_cellid = names(named_membership)[named_membership != '' & !is.na(named_membership)]
  named_membership = named_membership[assigned_cellid]
  
  test_stats = TabMC[TabMC$gamma == as.integer(tgamma),]
  
  obj_metacell = mcRigor_buildmc(obj_singlecell[, assigned_cellid], sc_membership = named_membership,
                                 fields = fields, aggregate_method = aggregate_method,
                                 add_testres = T, test_stats = test_stats, Thre = Thre)
  print(table(obj_metacell$mcRigor))
  
  return(list(obj_metacell = obj_metacell, 
              thre = Thre, 
              TabMC = TabMC,
              test_plot = test_res$test_plot,
              purity_plot = test_res$purity_plot))
}




#' @title mcRigor_OPTIMIZE
#'
#' @description Main functionality 2: To select the optimal hyperparameter for metacell partitioning 
#'
#'
#' @param obj_singlecell Seurat object of the original single cells.
#' @param cell_membership A dataframe, each column of which contains the metacell membership of single cells under a given gamma. 
#' The column names should be the corresponding gamma's. The row names should be the single cell id's.
#' @param Nrep The number of permutation repetitions we use when deriving the null.
#' @param assay_type The type of data assay yuo are using, depending on which different normalization would be used.
#' @param Gammas The candidate pool of granularity levels to consider in optimization
#' @param aggregate_method The method to aggregate single cell profiles into metacell profiles
#' @param output_file The directory for output saving
#' @param cor_method The method for gene correlation calculation description
#' @param gene_filter A proportion. Genes expressed lower than this proportion will be filtered out.
#' @param feature_use The number of genes to use in metacell testing.
#' @param prePro A boolean indicating whether to normalize obj_singlecell for preprocessing.
#' @param test_cutoff The test size for dubious metacell detection testing
#' @param thre_smooth A boolean indicating whether to smooth the threshold function
#' @param thre_bw If thre_smooth is True, thre_bw specifies the bandwidth for smoothing.
#' @param check_purity A boolean indicating whether to calculate the metacell purity of specific fields or not. 
#' @param fields A vector of the fields of interest, ex. celltype. It should be a subset of obj_singlecell's meta.data.
#' @param D_bw A boolean indicating whether to smooth the dubious rate with respect to metacell size
#' @param optim_method The method used for granularity level optimization. Default is trading off between sparsity and dubious rate
#' @param weight The weight for dubious rate in the tradeoff.
#' @param dub_rate If tradeoff is not used for optimization, what is highest acceptable dubious rate
#' @param draw  A boolean indicating whether to visualize the mcRigor results
#' @param pur_metric Can be NULL or a metadata variable name, ex. cell type.
#' @param step_save A boolean indicating whether to save the outputs step by step
#' 
#' 
#'
#' @return A list containing the following fields: 
#' \item{best_graining_level}{The optimal graining level selected by mcRigor}
#' \item{opt_metacell}{The metacell object build under the optimal graining level}
#' \item{scores}{A data frame containing the evaluation scores for each gamma}
#' \item{optim_plot}{The line plot to visualize the tradeoff for hyperparameter opimization.}
#' \item{thre}{The thresholds for dubious metacell detection}
#' \item{TabMC}{A dataframe containing the permutation results, elements to calculate the test statistics mcDiv and mcDiv null}
#'


#
mcRigor_OPTIMIZE <- function(obj_singlecell, 
                             cell_membership = NULL, 
                             assay_type = c('RNA', 'ATAC'),
                             Gammas = NULL,
                             aggregate_method = c('mean', 'sum', 'geom'),
                             output_file = NULL,
                             Nrep = 1,
                             gene_filter = 0.1, 
                             feature_use = 2000,
                             cor_method = c('pearson', 'spearman'),
                             prePro = T,
                             test_cutoff = 0.01, 
                             thre_smooth = T, 
                             thre_bw = 1/6,
                             D_bw = 10,
                             optim_method = c('tradeoff', 'dub_rate_large', 'dub_rate_small'),
                             weight = 0.5,
                             dub_rate=0.1, 
                             draw = T, 
                             pur_metric = NULL,
                             check_purity = T,
                             fields = NULL,
                             step_save = T){
  
  assay_type = match.arg(assay_type)
  cor_method = match.arg(cor_method)
  optim_method = match.arg(optim_method)
  aggregate_method = match.arg(aggregate_method)
  
  if (prePro){
    cat("Normalizing data...\n")
    if (assay_type == 'ATAC') obj_singlecell = Signac::RunTFIDF(obj_singlecell, verbose = F)
    else obj_singlecell = Seurat::NormalizeData(obj_singlecell, verbose = F)
  }
  
  if (is.null(fields)) {
    fields = sapply(colnames(obj_singlecell@meta.data),
                    function(X) {is.character(obj_singlecell[[X]][,1]) | is.factor(obj_singlecell[[X]][,1])})
    fields = names(fields[fields])
  }
  
  if (assay_type == 'ATAC'){
    feature_cutoff = paste0('q', floor((1 - feature_use / dim(obj_singlecell)[1]) * 100))
    obj_singlecell = Signac::FindTopFeatures(obj_singlecell, min.cutoff = feature_cutoff, verbose = F)
    top_gene = Seurat::VariableFeatures(obj_singlecell)
  } else {
    obj_singlecell = Seurat::FindVariableFeatures(obj_singlecell, nfeatures = feature_use, verbose = F)
    top_gene = Seurat::VariableFeatures(obj_singlecell)
  }
  
  if (is.null(Gammas)) Gammas = colnames(cell_membership)
  
  TabMC_list = vector(mode = 'list', length = length(Gammas))
  names(TabMC_list) = Gammas
  
  for (gamma in Gammas) {
    
    named_membership = cell_membership[[gamma]]
    names(named_membership) = rownames(cell_membership)
    assigned_cellid = names(named_membership)[named_membership != '' & !is.na(named_membership)]
    named_membership = named_membership[assigned_cellid]
    
    obj_metacell = mcRigor_buildmc(obj_singlecell[, assigned_cellid], sc_membership = named_membership,
                                   assay_type = assay_type, fields = fields, aggregate_method = 'sum',
                                   doNorm = F)
    obj_metacell[['ZeroRate']] = 1 - obj_metacell[[paste0('nFeature_', Seurat::DefaultAssay(obj_metacell))]] / nrow(obj_metacell)
    
    ##########
    
    metacell_id = colnames(obj_metacell)
    
    mc_stats = data.frame(matrix(NA, nrow = length(metacell_id), ncol = 3+2*Nrep))
    rownames(mc_stats) = metacell_id
    colnames(mc_stats) = c('size', 'T_org', 'T_colperm',
                           paste0('T_rowperm', 1:Nrep), paste0('T_bothperm', 1:Nrep))
    #
    cat('gamma =', gamma, '\n')
    pb <- utils::txtProgressBar(style=3)
    pbid = 0
    
    for (select.label in metacell_id) {
      
      select_sc = names(named_membership)[named_membership == select.label]
      mc_stats[select.label, 'size'] = length(select_sc)
      
      if (length(select_sc) == 1) next
      
      if(dim(obj_singlecell)[1] > dim(obj_singlecell)[2]){
        counts = Seurat::GetAssayData(object = obj_singlecell, layer = 'data')[top_gene, select_sc]
      } else {
        select_mc = obj_singlecell[, select_sc]
        counts = Seurat::GetAssayData(select_mc, layer = 'data')
      }
      
      counts = as.matrix(counts)
      # counts = apply(counts, 2, function(x) x/sqrt(sum(x^2)))
      
      counts_var = counts[apply(counts, 1, 
                                function(x) length(which(x>0)) > dim(counts)[2] * gene_filter), ]
      
      if (cor_method == 'pearson') dat = scale_cpp(t(counts_var))
      else {
        dat = apply(t(counts_var), 2, rank)
        dat = scale_cpp(dat)
      }
      
      p = dim(dat)[2]
      n = dim(dat)[1]
      
      mc_stats[select.label, 'T_org'] = mc_indpd_stats_cpp(dat)
      mc_stats[select.label, 'T_colperm'] = mc_indpd_stats_cpp(colwise_perm_cpp(dat))
      
      T_perm = sapply(1:Nrep, function(i) {
        x = rowwise_perm_cpp(dat)
        c(mc_indpd_stats_cpp(x), mc_indpd_stats_cpp(colwise_perm_cpp(x)))
      })
      mc_stats[select.label, paste0('T_rowperm', 1:Nrep)] = T_perm[1,]
      mc_stats[select.label, paste0('T_bothperm', 1:Nrep)] = T_perm[2,]
      
      pbid = pbid + 1
      utils::setTxtProgressBar(pb, pbid / length(metacell_id))
      
    }
    
    close(pb)
    
    purity_fields = colnames(obj_metacell@meta.data)[grep('_purity', colnames(obj_metacell@meta.data))]
    TabMC_list[[gamma]] =  cbind(as.numeric(gamma), 
                                 obj_metacell[['ZeroRate']],
                                 obj_metacell[[purity_fields]], 
                                 mc_stats)
    rownames(TabMC_list[[gamma]]) = rownames(mc_stats)
    
    if (step_save){
      TabMC = do.call(rbind, unname(TabMC_list))
      colnames(TabMC)[1] = 'gamma'
      colnames(TabMC)[2] = 'ZeroRate'
      
      TabMC$TT_div = TabMC$T_org / TabMC$T_colperm
      TabMC$TT_div[is.na(TabMC$TT_div)] = 1
      
      if (is.null(output_file)) 
        output_file = paste0('Tabmc_', assay_type, '_optimize.RData')
      save(TabMC, file = output_file)
    }
    
  }
  
  TabMC = do.call(rbind, unname(TabMC_list))
  colnames(TabMC)[1] = 'gamma'
  colnames(TabMC)[2] = 'ZeroRate'
  
  TabMC$TT_div = TabMC$T_org / TabMC$T_colperm
  TabMC$TT_div[is.na(TabMC$TT_div)] = 1
  
  if (is.null(output_file)) 
    output_file = paste0('Tabmc_', assay_type, '_optimize.RData')
  save(TabMC, file = output_file)
  
  ###########
  
  test_res = mcRigor_threshold(TabMC, 
                               test_cutoff = test_cutoff,
                               thre_smooth = thre_smooth, thre_bw = thre_bw,
                               draw = draw, pur_metric = pur_metric)
  Thre = test_res$threshold
  TabMC = test_res$TabMC
  
  if (is.null(output_file)) 
    output_file = paste0('Tabmc_', assay_type, '_optimize.RData')
  save(TabMC, file = output_file)
  
  tradeoff_res = mcRigor_tradeoff(TabMC, 
                                  threshold = Thre,
                                  optim_method = optim_method, 
                                  D_bw = D_bw,
                                  weight = weight,
                                  dub_rate = dub_rate, 
                                  draw = draw)
  opt_gamma = tradeoff_res$optimized$gamma
  
  ##
  
  opt_named_membership = cell_membership[[as.character(opt_gamma)]]
  names(opt_named_membership) = rownames(cell_membership)
  assigned_cellid = names(opt_named_membership)[opt_named_membership != '' & !is.na(opt_named_membership)]
  opt_named_membership = opt_named_membership[assigned_cellid]
  
  test_stats = TabMC[TabMC$gamma == opt_gamma,]
  
  obj_metacell = mcRigor_buildmc(obj_singlecell[, assigned_cellid], sc_membership = opt_named_membership,
                                 fields = fields, aggregate_method = aggregate_method,
                                 add_testres = T, test_stats = test_stats, Thre = Thre)
  cat('optimal gamma =', opt_gamma, '\n')
  print(table(obj_metacell$mcRigor))
  
  return(list(best_graining_level = opt_gamma, 
              opt_metacell = obj_metacell, 
              scores = tradeoff_res$scores,
              thre = Thre, 
              TabMC = TabMC,
              optim_plot = tradeoff_res$optim_plot))
}






#' @title mcRigor_threshold
#'
#' @description A building block of the main functions. To derive the thresholds for detecting dubious metacells based on the output permutation results (TabMC)
#'
#'
#' @param TabMC A dataframe containing the permutation results. Saved in the previous steps
#' @param test_cutoff The test size for dubious metacell detection testing
#' @param thre_smooth A boolean indicating whether to smooth the threshold function
#' @param thre_bw If thre_smooth is True, what is the bandwidth for smoothing
#' @param draw  A boolean indicating whether to visualize the mcRigor results
#' @param palpha Point alpha value for transparency in drawing.
#' @param org_color The colors indicating metacell purities or other interested factors
#' @param null_color The color for the null.
#' @param pur_metric Name of the covariate that we want to compute purity on. Can be NULL or a metadata variable name, ex. cell type.
#' 
#'
#' @return A list containing the following fields: 
#' \item{threshold}{The thresholds for dubious metacell detection}
#' \item{TabMC}{A dataframe containing the permutation results and the testing results given by mcRigor}
#' \item{test_plot}{The scatter plots demonstrating the mcDiv values and the obtained thresholds for dubious metacell detection}
#' \item{purity_plot}{A violin plot showing the purity distribution of the pur_metric covariate in dubious metacells and trustworthy metacells}
#'


mcRigor_threshold <- function(TabMC, 
                              test_cutoff = 0.01, 
                              thre_smooth = T, 
                              thre_bw = 1/6,
                              draw = T,
                              palpha = 1,
                              org_color = c('red', 'orange','yellow', 'lightblue'),  # c('#CCA453')
                              null_color = '#666666',
                              pur_metric = NULL) {
  
  TabMC$depctr_TT_div = TabMC$T_rowperm1 / TabMC$T_bothperm1
  TabMC$depctr_TT_div[is.na(TabMC$depctr_TT_div)] = 1
  
  rowperm_col = grep('T_rowperm', colnames(TabMC))
  bothperm_col = grep('T_bothperm', colnames(TabMC))
  # TabMC$ctrl = apply(TabMC, 1, function(x) max(x[rowperm_col] / x[bothperm_col]))
  
  if (draw) {
    
    ggplot2::theme_set(ggplot2::theme_classic())
    
    if (!is.null(pur_metric) && !pur_metric %in% colnames(TabMC)) pur_metric = colnames(TabMC)[grep(pur_metric, colnames(TabMC))[1]]
    if (is.null(pur_metric)) pur_metric = colnames(TabMC)[grep('type', colnames(TabMC))[1]]
    if (is.na(pur_metric))  pur_metric = colnames(TabMC)[grep('_purity', colnames(TabMC))[length(grep('_purity', colnames(TabMC)))]]
    if (is.na(pur_metric)) {
      pur_metric = 'ident'
      TabMC[[pur_metric]] = 'all'
    }
    
    pur_min = min(TabMC[[pur_metric]])
    
    ymin = min(min(TabMC$TT_div), min(TabMC$depctr_TT_div))
    ymax = max(max(TabMC$TT_div), max(TabMC$depctr_TT_div))
    
    p1 = ggplot2::ggplot(TabMC, ggplot2::aes_string(x='size', y='depctr_TT_div')) +
      ggplot2::geom_point(mapping = ggplot2::aes_string(color=pur_metric), alpha=palpha, col = null_color) + 
      # ggplot2::scale_color_gradientn(name=pur_metric, limits = c(pur_min, 1), colors = org_color) +
      ggplot2::theme(legend.position = 'none')  +
      ggplot2::ylab('mcDiv_null') + ggplot2::xlab('metacell size') +
      ggplot2::ylim(c(ymin, ymax)) +
      # ggplot2::ylim(c(0.95, 1.6)) +
      ggplot2::xlim(c(0, min(max(TabMC$size), 500)))
    
    p2 = ggplot2::ggplot(TabMC, ggplot2::aes_string(x='size', y='TT_div')) +
      ggplot2::geom_point(mapping = ggplot2::aes_string(color=pur_metric), alpha=palpha) + 
      # scale_color_viridis_c(name = 'purity', option = 'C') +
      ggplot2::scale_color_gradientn(limits = c(pur_min, 1), colors = org_color) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), 
                     legend.key.width = ggplot2::unit(0.3,'cm'), legend.text = ggplot2::element_text(size = 7)) +
      ggplot2::ylab('mcDiv') + ggplot2::xlab('metacell size') +
      ggplot2::ylim(c(ymin, ymax)) +
      # ggplot2::ylim(c(0.95, 1.6)) +
      ggplot2::xlim(c(0, min(max(TabMC$size), 500)))
    suppressWarnings(p_legend <- cowplot::get_legend(p2))
    p2 = p2 + ggplot2::theme(legend.position = 'none')
    
    p12 = p2 + 
      ggplot2::geom_point(mapping = ggplot2::aes_string(x='size', y='depctr_TT_div', color=pur_metric), alpha=0.1, col = null_color) +
      ggplot2::xlab('metacell size') + ggplot2::ylab('statistics')
    
    # p2+p1
  }
  
  
  Thre = as.data.frame(NULL)
  
  for (size in unique(TabMC$size)) {
    
    if (size == 1) next
    
    mcs = TabMC[TabMC$size==size,]
    
    thre = stats::quantile(unlist(mcs[, rowperm_col] / mcs[, bothperm_col]), probs = 1-test_cutoff)
    
    Thre = rbind(Thre, c(size, thre))
  }
  
  colnames(Thre) = c('size', 'thre')
  Thre = Thre[order(Thre$size),]
  
  if (thre_smooth) {
    a=stats::lowess(x=Thre$size, y=Thre$thre, f=thre_bw)
    # plot(Thre$size, Thre$thre)
    # points(a$x, a$y)
    Thre$thre = a$y
    Thre$size = a$x
  }
  
  if (draw) {
    pthre1 = p1 + ggplot2::geom_path(data = Thre, mapping = ggplot2::aes_string(x='size', y='thre'), color='red') 
    pthre2 = p2 + ggplot2::geom_path(data = Thre, mapping = ggplot2::aes_string(x='size', y='thre'), color='red') 
    # pthre2+pthre1
    pthre12 = p12 + ggplot2::geom_path(data = Thre, mapping = ggplot2::aes_string(x='size', y='thre'), color='red') 
  }
  
  TabMC$mcRigor = 'trustworthy'
  for (mcid in 1:dim(TabMC)[1]) {
    if (TabMC$size[mcid] >1 && TabMC$TT_div[mcid] > Thre$thre[Thre$size == TabMC$size[mcid]]) TabMC$mcRigor[mcid] = 'dubious'
    # if (TabMC$size[mcid] >1 && all(TabMC$TT_div[mcid] > TabMC[mcid, rowperm_col] / TabMC[mcid, bothperm_col]) ) TabMC$mcRigor[mcid] = 'dubious'
    
  }
  # print(table(TabMC$mcRigor))
  
  pvioin = NULL
  if (draw) {
    pviolin = ggplot2::ggplot(TabMC) + ggplot2::geom_violin(mapping = ggplot2::aes_string(y=pur_metric, x='mcRigor', fill='mcRigor')) +
      ggplot2::theme(legend.position = 'none') + ggplot2::xlab(NULL) + ggplot2::ylab(pur_metric)
    # pviolin = ggplot2::ggplot(TabMC) + geom_boxplot(mapping = ggplot2::aes_string(y=pur_metric, x='mcRigor', fill='mcRigor')) +
    #   ggplot2::theme(legend.position = 'none') + ggplot2::xlab(NULL) + ggplot2::ylab(pur_metric) +
    #   ggplot2::scale_fill_manual(values = c('#310690', '#F0F921'))
  }
  
  ###
  
  test_plot = NULL
  if (draw) {
    test_plot = cowplot::plot_grid(cowplot::plot_grid(pthre2, pthre1, nrow = 1),
                                   cowplot::plot_grid(p_legend, nrow = 1), 
                                   nrow = 1, rel_widths = c(1,0.1))
  }
  
  return(list(threshold = Thre, 
              TabMC = TabMC,
              test_plot = test_plot,
              purity_plot = pviolin))
}





#' @title mcRigor_tradeoff
#'
#' @description A building block of the main functions. To evaluate each metacell partition and optimize metacell partitioning based on the output permutation results (TabMC) and thresholds (threshold)
#'
#'
#' @param TabMC A dataframe containing the permutation results. Saved in the previous steps
#' @param threshold A dataframe containing the dubious metacell detection thresholds given by mcRigor_threshold
#' @param D_bw A boolean indicating whether to smooth the dubious rate with respect to metacell size
#' @param optim_method The method used for granularity level optimization. Default is trading off between sparsity and dubious rate
#' @param weight The weight for dubious rate in the tradeoff.
#' @param dub_rate If tradeoff is not used for optimization, what is highest acceptable dubious rate
#' @param draw  A boolean indicating whether to visualize the mcRigor results
#' 
#'
#' @return A list containing the following fields: 
#' \item{optimized}{The optimization results, containing the optimal gamma and its corresponding Sore}
#' \item{scores}{A data frame containing the evaluation scores for each gamma}
#' \item{optim_plot}{The line plot to visualize the tradeoff for hyperparameter opimization.}
#'


mcRigor_tradeoff <- function(TabMC, 
                             threshold,
                             D_bw = 10,
                             optim_method = c('tradeoff', 'dub_rate_large', 'dub_rate_small'),
                             dub_rate=0.1, 
                             weight = 0.5,
                             draw = T) {
  
  optim_method = match.arg(optim_method)
  
  Thre = threshold
  
  DD <- data.frame(gamma=sort(unique(TabMC$gamma)))
  # DD <- data.frame(gamma=DD$gamma[DD$gamma %% 5 == 0])
  DD$D <- 0
  DD$ZeroRate <- 0
  
  for (gamma in DD$gamma) {
    
    mcs = TabMC[TabMC$gamma==gamma & TabMC$mcRigor == 'trustworthy',]
    dub_mc = TabMC[TabMC$gamma==gamma & TabMC$mcRigor == 'dubious',]
    
    DD$D[DD$gamma == gamma] = sum(dub_mc$size) / (sum(mcs$size) + sum(dub_mc$size)) 
    
    DD$ZeroRate[DD$gamma == gamma] = mean(TabMC[TabMC$gamma == gamma, 'ZeroRate'])
    
  }
  
  DD = DD[!is.nan(DD$D),]
  if (max(DD$gamma) - min(DD$gamma) > 5){
    DD$D = stats::lowess(x = DD$gamma, y = DD$D, f = D_bw / (max(DD$gamma) - min(DD$gamma)))$y
  }
  # DD$score = (DD$ZeroRate * (max(DD$D) - min(DD$D)) + DD$D * (max(DD$ZeroRate) - min(DD$ZeroRate))) /
  #   (max(DD$D) - min(DD$D) + max(DD$ZeroRate) - min(DD$ZeroRate))
  if (is.null(weight)) {
    weight = (max(DD$ZeroRate) - min(DD$ZeroRate)) / (max(DD$D) - min(DD$D) + max(DD$ZeroRate) - min(DD$ZeroRate))
  }
  DD$score = DD$ZeroRate * (1-weight) + DD$D * weight
  # DD$score = lowess(x = DD$gamma, y = DD$score, f = 7 / (max(DD$gamma) - min(DD$gamma)))$y  # bandwidth around 5-10, f = 1/20 or f = 10 / (max(DD$gamma) - min(DD$gamma)) 
  DD$score = 1 - DD$score
  
  opt_gamma = switch(optim_method,
                     tradeoff = DD$gamma[which.max(DD$score)],
                     dub_rate_large = max(DD$gamma[DD$D < dub_rate]),
                     dub_rate_small = DD$gamma[which(DD$D >= dub_rate)[1] - 1])
  
  optim = data.frame(gamma = opt_gamma, D = DD$D[DD$gamma==opt_gamma], 
                     score = DD$score[DD$gamma == opt_gamma], 
                     cross_score = 1 - 1/2 * DD$D[DD$gamma==opt_gamma] - 1/2 * DD$ZeroRate[DD$gamma==opt_gamma])
  
  optimized = as.list(optim)
  
  optim_plot = NULL
  
  if (draw) {
    
    if (optim_method == 'tradeoff') {
      pelbow = ggplot2::ggplot(DD, ggplot2::aes(x=gamma, y=D)) + ggplot2::geom_path(ggplot2::aes(col = 'D')) +
        # ggplot2::geom_point(data = optim, mapping = ggplot2::aes_string(x='gamma', y='D'), col='red', size=3) +
        ggplot2::annotate('text', x=optim$gamma, y=-0.03, label=optim$gamma, color='red', size = 3) +
        ggplot2::geom_path(ggplot2::aes(x = gamma, y = ZeroRate, col = 'ZeroRate')) + 
        ggplot2::geom_path(ggplot2::aes(x = gamma, y = score, col = 'Score')) +
        ggplot2::scale_color_manual(name = NULL, values = c('Score' = 'darkred', 'D' = 'darkblue', 'ZeroRate' = 'darkgreen')) +
        ggplot2::geom_vline(xintercept = optim$gamma, col = 'red', linetype = 'dashed') +
        ggplot2::ylab(' ') +
        ggplot2::coord_cartesian(ylim = c(0,1), clip = 'off') +
        ggplot2::theme_light()+ ggplot2::theme(aspect.ratio = 2/3)
    } else {
      pelbow = ggplot2::ggplot(DD, ggplot2::aes_string(x='gamma', y='D')) + ggplot2::geom_path() +
        ggplot2::geom_point(data = optim, mapping = ggplot2::aes_string(x='gamma', y='D'), col='red', size=3) +
        ggplot2::annotate('text', x=optim$gamma, y=optim$D-0.05, label=optim$gamma, color='red', size = 3) +
        ggplot2::coord_cartesian(ylim = c(0,1), clip = 'off') +
        ggplot2::theme_light()+ ggplot2::theme(aspect.ratio = 2/3)
    }
    
    
    optim_plot = pelbow
    
  }
  
  return(list(optimized = optimized,
              scores = DD,
              optim_plot = optim_plot))
}







#' @title mcRigor_buildmc
#'
#' @description To build an metacell object from the given metacell partitioning (membership)
#'
#'
#' @param obj_singlecell the Seurat object of single cells
#' @param sc_membership A named vector (or dataframe) of the metacell membership of single cells
#' @param aggregate_method The method to aggregate single cell profiles into metacell profiles
#' @param assay_type The type of data assay yuo are using, depending on which different normalization would be used.
#' @param doAssign A bool indicating whether to assign covariates to metacells or not
#' @param fields A vector of covariate names to assign
#' @param covariate_method Method to define the most abundant cell covariate within metacells. Available: "jaccard", "relative", "absolute" (default).
#' \itemize{
#'   \item jaccard - assign metacell to covariate with the maximum jaccard coefficient (recommended)
#'   \item relative - assign metacell to covariate with the maximum relative abundance (normalized by cluster size), may result in assignment of metacells to poorly represented (small) covariate due to normalization
#'   \item absolute - assign metacell to covariate with the maximum absolute abundance within metacell, may result in disappearance of poorly represented (small) clusters
#' }
#' @param purity_method method to compute metacell purity.
#' \code{"max_proportion"} if the purity is defined as a proportion of the most abundant covariate (cell type) within super-cell or
#' \code{"entropy"} if the purity is defined as the Shanon entropy of the covariates metacell consists of.
#' @param doNorm A bool indicating whether to perform normalization for the metacell object (obj_metacell)
#' 
#' @param add_testres A bool indicating whether to add the mcRigor results (dubious or trustworthy) as part of obj_metacell's metadata
#' @param test_stats If add_testres = True, this argument is needed. Usually should be TabMC, an output from previous steps.
#' @param Thre The threshold for dubious metacell detection. If not inputed, it will be computed based on test_stats.
#' @param test_cutoff The test size for dubious metacell detection testing
#' @param cor_method The method for gene correlation calculation description
#' @param gene_filter A proportion. Genes expressed lower than this proportion will be filtered out.
#' @param feature_use The number of genes to use in metacell testing.
#' @param prePro A boolean indicating whether to normalize obj_singlecell for preprocessing.

#' 
#' @return a Seurat object of the metacells
#'
#' @export
#'

mcRigor_buildmc <- function(obj_singlecell, 
                            sc_membership = NULL, 
                            assay_type = c('RNA', 'ATAC'), 
                            doNorm = T,
                            aggregate_method = c('mean', 'sum', 'geom'),
                            doAssign = T, 
                            fields = NULL,
                            covariate_method = 'absolute', 
                            purity_method = 'max_proportion',
                            add_testres = F, test_stats = NULL, Thre = NULL, test_cutoff = 0.01,
                            prePro = F, feature_use = 2000, gene_filter = 0.1, 
                            cor_method = c('pearson', 'spearman')) {
  
  if (is.null(sc_membership)) stop('Please provide the metacell membership of single cells.')
  if (is.null(names(sc_membership))) warning('The metacell memberships are not matched to the single cell ids.')
  
  assay_type = match.arg(assay_type)
  aggregate_method = match.arg(aggregate_method)
  
  obj_singlecell[['Metacell']] = sc_membership
  sc_membership = obj_singlecell$Metacell
  
  if (aggregate_method == 'sum') {
    counts_metacell = Seurat::AggregateExpression(obj_singlecell,
                                                  return.seurat = F,
                                                  group.by = 'Metacell',
                                                  verbose  = F)[[1]]
    if (is.numeric(sc_membership[1])) {
      mc_names = sapply(colnames(counts_metacell), function(x) substr(x, start = 2, stop = 100000))
      colnames(counts_metacell) = as.numeric(mc_names)
    }
  } else if (aggregate_method == 'mean') {
    counts_sc = Seurat::GetAssayData(obj_singlecell, layer = 'counts')
    counts_metacell = t(apply(counts_sc, 1, function(x) tapply(x, sc_membership, mean)))
  } else {
    counts_sc = Seurat::GetAssayData(obj_singlecell, layer = 'counts')
    counts_metacell = t(apply(counts_sc, 1, function(x) tapply(x, sc_membership, function(y) {
      exp(mean(log(1 + y))) - 1
    })))
  }
  
  if (assay_type == 'ATAC') {
    
    suppressWarnings(obj_metacell <- Signac::CreateChromatinAssay(counts = counts_metacell,
                                                                  sep = c('-', '-')))
    suppressWarnings(obj_metacell <- Seurat::CreateSeuratObject(counts = obj_metacell,
                                                                assay = 'peaks'))
    # obj_metacell[['peaks']]@meta.features = obj_singlecell[['peaks']]@meta.features
    if (doNorm) obj_metacell = Signac::RunTFIDF(obj_metacell, verbose = F)
    
  } else {
    
    suppressWarnings(obj_metacell <- Seurat::CreateSeuratObject(counts = counts_metacell))
    if (doNorm) obj_metacell = Seurat::NormalizeData(obj_metacell, verbose = F)
    
  }
  
  
  if (doAssign) {
    
    if (is.null(fields)) {
      fields = sapply(colnames(obj_singlecell@meta.data),
                      function(X) {is.character(obj_singlecell[[X]][,1]) | is.factor(obj_singlecell[[X]][,1])})
      fields = names(fields[fields])
      fields = fields[sapply(fields, function(X) length(table(obj_singlecell[[X]])) > 1)]
    }
    
    # cat("Assign metadata to metacells and compute purities...\n")
    for (f in fields) {
      
      if (!all(is.na(obj_singlecell[[f]][,1]))) {
        obj_metacell[[f]] = mcRigor_covariate(sc_covariate = obj_singlecell[[f]],
                                              sc_membership = obj_singlecell[['Metacell']],
                                              method = covariate_method)
        if(length(unique(obj_singlecell[[f]][,1]))> 1) {
          obj_metacell[[paste0(f,"_purity")]] = mcRigor_purity(sc_covariate = obj_singlecell[[f]],
                                                               sc_membership = obj_singlecell[['Metacell']],
                                                               method = purity_method)
        }
      }
    }
  }
  
  cl.gr = table(sc_membership)
  obj_metacell[['size']] = c(cl.gr)
  
  obj_metacell@misc$cell_membership = obj_singlecell[['Metacell']]
  
  if ('umap' %in% Seurat::Reductions(obj_singlecell) || 'tsne' %in% Seurat::Reductions(obj_singlecell)) {
    reduc_key = ifelse('umap' %in% Seurat::Reductions(obj_singlecell), 'umap', 'tsne')
    sc_embed = Seurat::Embeddings(obj_singlecell, reduction = reduc_key)
    mc_embed = stats::aggregate(sc_embed, by = list(sc_membership), FUN = mean)
    rownames(mc_embed) = mc_embed[[1]]
    obj_metacell[[names(mc_embed)[2]]] = mc_embed[2]
    obj_metacell[[names(mc_embed)[3]]] = mc_embed[3]
  }
  
  if (add_testres) {
    if (is.null(test_stats)) 
      stop('To add test results to the metacell object, please provide the calculated test statistics (test_stats).')
    if (is.null(Thre)) {
      Thre = mcRigor_threshold(test_stats, test_cutoff = test_cutoff)$threshold
      cat('Testing thresholds (Thre) not provided. Calculating thresholds with defaults...')
    }
    
    if (!all(colnames(obj_metacell) %in% rownames(test_stats))){
      obj_metacell$TT_div = NA
      if (prePro){
        cat("Normalizing sc data...\n")
        if (assay_type == 'ATAC') obj_singlecell = Signac::RunTFIDF(obj_singlecell, verbose = F)
        else obj_singlecell = Seurat::NormalizeData(obj_singlecell, verbose = F)
      }
      
      if (assay_type == 'ATAC'){
        feature_cutoff = paste0('q', floor((1 - feature_use / dim(obj_singlecell)[1]) * 100))
        obj_singlecell = Signac::FindTopFeatures(obj_singlecell, min.cutoff = feature_cutoff, verbose = F)
        top_gene = Seurat::VariableFeatures(obj_singlecell)
      } else {
        obj_singlecell = Seurat::FindVariableFeatures(obj_singlecell, nfeatures = feature_use, verbose = F)
        top_gene = Seurat::VariableFeatures(obj_singlecell)
      }
    }
    
    obj_metacell$mcRigor = 'trustworthy'
    obj_metacell@misc$cell_membership$mcRigor_sc = 'trustworthy'
    
    for (mcid in colnames(obj_metacell)) {
      
      if (mcid %in% rownames(test_stats)){
        
        if (obj_metacell$size[mcid] > 1 && 
            test_stats[mcid, 'TT_div'] > Thre$thre[Thre$size == obj_metacell$size[mcid]]) {
          obj_metacell$mcRigor[mcid] = 'dubious'
          obj_metacell@misc$cell_membership$mcRigor_sc[obj_metacell@misc$cell_membership$Metacell == mcid] = 'dubious'
        }
        
      } else {
        
        select_sc = names(sc_membership)[sc_membership == mcid]
        
        if (length(select_sc) == 1) next
        
        if(dim(obj_singlecell)[1] > dim(obj_singlecell)[2]){
          counts = Seurat::GetAssayData(object = obj_singlecell, layer = 'data')[top_gene, select_sc]
        } else {
          select_mc = obj_singlecell[, select_sc]
          counts = Seurat::GetAssayData(select_mc, layer = 'data')
        }
        
        counts = as.matrix(counts)
        # counts = apply(counts, 2, function(x) x/sqrt(sum(x^2)))
        
        counts_var = counts[apply(counts, 1, 
                                  function(x) length(which(x>0)) > dim(counts)[2] * gene_filter), ]
        
        cor_method = match.arg(cor_method)
        if (cor_method == 'pearson') dat = scale_cpp(t(counts_var))
        else {
          dat = apply(t(counts_var), 2, rank)
          dat = scale_cpp(dat)
        }
        
        p = dim(dat)[2]
        n = dim(dat)[1]
        
        mc_T_org = mc_indpd_stats_cpp(dat)
        mc_T_colperm = mc_indpd_stats_cpp(colwise_perm_cpp(dat))
        mc_TT_div = mc_T_org / mc_T_colperm
        
        obj_metacell$TT_div[mcid] = mc_TT_div
        
        if (obj_metacell$size[mcid] %in% Thre$size) {
          if (mc_TT_div > Thre$thre[Thre$size == obj_metacell$size[mcid]]) {
            obj_metacell$mcRigor[mcid] = 'dubious'
            obj_metacell@misc$cell_membership$mcRigor_sc[obj_metacell@misc$cell_membership$Metacell == mcid] = 'dubious'
          }
        } else {
          temp_thre = stats::approx(x = Thre$size, y = Thre$thre, xout = obj_metacell$size[mcid], 
                                    yleft = min(Thre$thre), yright = max(Thre$thre))$y
          if (mc_TT_div > temp_thre) {
            obj_metacell$mcRigor[mcid] = 'dubious'
            obj_metacell@misc$cell_membership$mcRigor_sc[obj_metacell@misc$cell_membership$Metacell == mcid] = 'dubious'
          }
        }
        
      }
    }
    # print(table(obj_metacell$mcRigor))
    
  }
  
  return(obj_metacell)
  
}








#' @title mcRigor_covariate
#'
#' @description To assist metacell object building in determining the metacell covariates (metadata)
#'
#'
#' @param sc_covariate a vector (or dataframe) of single cell covariates
#' @param sc_membership a vector (or dataframe) of metacell memberships of single cells 
#' @param method method to define the most abundant cell covariate within metacells. Available: "jaccard" (default), "relative", "absolute".
#' \itemize{
#'   \item jaccard - assign metacell to covariate with the maximum jaccard coefficient (recommended)
#'   \item relative - assign metacell to covariate with the maximum relative abundance (normalized by cluster size), may result in assignment of metacells to poorly represented (small) covariate due to normalization
#'   \item absolute - assign metacell to covariate with the maximum absolute abundance within metacell, may result in disappearance of poorly represented (small) clusters
#' }
#'
#' @return a vector of assigned metacell covariates
#'
#' @export
#'


mcRigor_covariate <- function(sc_covariate, sc_membership, 
                              method = c("jaccard", "relative", "absolute")){
  
  if (is.list(sc_membership) && is.list(sc_covariate)) {
    sc_membership = sc_membership[rownames(sc_covariate), 1]
    sc_covariate = sc_covariate[,1]
  } else if (is.list(sc_membership)) {
    sc_membership = sc_membership[, 1]
  } else if (is.list(sc_covariate)) {
    sc_covariate = sc_covariate[,1]
  }
  
  cl.gr            <- table(sc_covariate, sc_membership)
  covariate.size     <- as.numeric(table(sc_covariate))
  group.size       <- as.numeric(table(sc_membership))
  
  if(is.null(method[1]) | is.na(method[1]) | is.nan(method[1])){
    stop(paste("Please specify method: jaccard (recommended), relative or absolute"))
  }
  
  if(method[1] == "jaccard"){
    
    cl.gr          <- as.matrix(cl.gr)
    jaccard.mtrx   <- cl.gr
    
    for(i in rownames(cl.gr)){
      for(j in colnames(cl.gr)){
        jaccard.mtrx[i,j] <- cl.gr[i,j] / (sum(cl.gr[i,]) +  sum(cl.gr[,j]) - cl.gr[i,j])
      }
    }
    res              <- apply(jaccard.mtrx, 2, function(x){names(x)[which.max(x)]})
    
  } else if(method[1] == "relative"){
    cl.gr            <- sweep(cl.gr, 1, covariate.size, "/")
    res              <- apply(cl.gr, 2, function(x){names(x)[which.max(x)]})
    
  } else if(method[1] == "absolute"){
    res              <- apply(cl.gr, 2, function(x){names(x)[which.max(x)]})
    
  } else {
    stop(paste("Unknown value of method (", method[1] , ")", "please, use: jaccard, relative or absolute"))
  }
  
  
  return(res)
}






#' @title mcRigor_purity
#'
#' @description To assist metacell object building in compute the purity of metacells with respect to each covariate
#'
#'
#' @param sc_covariate vector of single cell covariates
#' @param sc_membership vector of assignment of single-cell data to metacells
#' @param method method to compute metacell purity.
#' \code{"max_proportion"} if the purity is defined as a proportion of the most abundant covariate (cell type) within super-cell or
#' \code{"entropy"} if the purity is defined as the Shanon entropy of the covariates metacell consists of.
#'
#' @return a vector of metacell purity, which is defined as:
#' - proportion of the most abundant covariate within metacell for \code{method = "max_proportion"} or
#' - Shanon entropy for \code{method = "entropy"}.
#' With 1 meaning that metacell consists of single cells from one covariate (reference assignment)
#'
#' @export
#'

mcRigor_purity <- function(sc_covariate, sc_membership,
                           method = c("max_proportion", "entropy")[1]){
  
  if(!(method %in% c("max_proportion", "entropy"))){
    stop(paste("Method", method, "is not known. The available methods are:", paste(method, collapse = ",")))
  }
  
  if (is.list(sc_membership) && is.list(sc_covariate)) {
    sc_membership = sc_membership[rownames(sc_covariate), 1]
    sc_covariate = sc_covariate[,1]
  } else if (is.list(sc_membership)) {
    sc_membership = sc_membership[, 1]
  } else if (is.list(sc_covariate)) {
    sc_covariate = sc_covariate[,1]
  }
  
  cl.gr            <- table(sc_covariate, sc_membership)
  
  switch(method,
         
         entropy = {
           res <- apply(cl.gr, 2, entropy::entropy)
         },
         
         max_proportion = {
           cluster.size     <- as.numeric(table(sc_covariate))
           group.size       <- as.numeric(table(sc_membership))
           
           Ng               <- length(group.size)
           # group.max.cl     <- rep(0, Ng)
           
           
           cl.gr            <- sweep(cl.gr, 2, group.size, "/")
           
           res              <- apply(cl.gr, 2, max)
         },
         
         {
           stop(paste("Method", method, "is not known. The available methods are:", paste(method, collapse = ",")))
         }
  )
  
  
  return(res)
}




#' @title mcRigor_projection
#'
#' @description To visualize a given metacell partitioning by projecting metacells onto the single cell embedding space.
#'
#'
#' @param obj_singlecell the Seurat object of single cells
#' @param sc_membership A named vector (or dataframe) of the metacell membership of single cells
#' @param sc.reduction The single cell reduction method to produce single cell embeddings
#' @param dims The dimensions to use in the single cell embeddings
#' @param metric The variable that determines the sizes of dots representing the metacells
#' @param color_field The variable based on which to color the cells.
#' 
#' @param add_testres A bool indicating whether to add the mcRigor results (dubious or trustworthy) as part of obj_metacell's metadata
#' @param test_stats If add_testres = True, this argument is needed. Usually should be TabMC, an output from previous steps.
#' @param Thre The threshold for dubious metacell detection. If not inputed, it will be computed based on test_stats.
#' @param test_cutoff The test size for dubious metacell detection testing
#'
#'
#'
#' @return a vector of metacell purity, which is defined as:
#' - proportion of the most abundant covariate within metacell for \code{method = "max_proportion"} or
#' - Shanon entropy for \code{method = "entropy"}.
#' With 1 meaning that metacell consists of single cells from one covariate (reference assignment)
#'
#' @export
#'

mcRigor_projection <- function(obj_singlecell, sc_membership = NULL,
                               sc.reduction = c("umap", "tsne", "pca"),
                               dims = c(1, 2),
                               metric = "size",
                               add_testres = F, test_stats = NULL, Thre = NULL, test_cutoff = 0.01,
                               color_field = NULL,
                               cpalette = c('#13678A', '#45C4B0', '#9AEBA3', '#BF847E', '#F2C12E', "#7E57C2", '#FACFCE', "#9E9D24",
                                            '#DAFDBA', '#86ABD4', "#42A5F5", "#546E7A", "#D4E157", "#76FF03", "#6D4C41", "#004D40",
                                            "#AB47BC", "#D81B60"),
                               sc.alpha = 0.5,
                               mc.alpha = 1,
                               pt_size = 1,
                               axis_lab = F,
                               max_mcsize = 1000,
                               continuous_metric = T,
                               dub_mc.label = F,
                               dub_mc_test.label = F,
                               label_text = F
) {
  
  if (is.null(sc_membership)) stop('Please provide the metacell membership of single cells.')
  if (is.null(names(sc_membership))) warning('The metacell memberships are not matched to the single cell ids.')
  
  sc.reduction = match.arg(sc.reduction)
  
  ggplot2::theme_set(ggplot2::theme_classic())
  
  obj_singlecell[['Metacell']] = sc_membership
  sc_membership = obj_singlecell$Metacell
  
  
  if(is.character(sc.reduction)){
    # if sc_reduction does not exist compute pca and run UMAP:
    if(is.null(obj_singlecell@reductions[[sc.reduction]])){
      cat("Low dimensionnal embessing not found in obj_singlecell")
      cat("Computing PCA ...")
      obj_singlecell <- Seurat::NormalizeData(obj_singlecell, verbose = F)
      obj_singlecell <- Seurat::FindVariableFeatures(obj_singlecell, verbose = F)
      obj_singlecell <- Seurat::ScaleData(obj_singlecell, verbose = F)
      obj_singlecell <- Seurat::RunPCA(obj_singlecell, verbose = F)
      cat("Running UMAP ...")
      obj_singlecell <- Seurat::RunUMAP(obj_singlecell, reduction = "pca", dims = c(1:20), verbose = F)
      scCoord <- Seurat::Embeddings(obj_singlecell@reductions[["umap"]])
    } else{
      scCoord <- Seurat::Embeddings(obj_singlecell@reductions[[sc.reduction]])
    } 
  } else{
    scCoord <- sc.reduction
  }  
  
  centroids <- stats::aggregate(scCoord~sc_membership, scCoord, mean) #should be taken from object slot
  
  add_testres = dub_mc_test.label
  obj_metacell = mcRigor_buildmc(obj_singlecell, sc_membership, doNorm = F,
                                 add_testres = add_testres, test_stats = test_stats, Thre = Thre, test_cutoff = test_cutoff)
  
  centroids = centroids[match(colnames(obj_metacell), centroids[,1]),]
  
  centroids[[metric]] <- obj_metacell[[metric]][,1]
  
  if ('celltype_purity' %in% colnames(obj_metacell@meta.data) & dub_mc.label) {
    centroids[['dub_mc']] = NA
    centroids[['dub_mc']][obj_metacell[['celltype_purity']]<0.75] = colnames(obj_metacell)[obj_metacell[['celltype_purity']]<0.75]
  }
  
  if ('metacell_purity' %in% colnames(obj_metacell@meta.data) & dub_mc.label) {
    centroids[['dub_mc']] = NA
    centroids[['dub_mc']][obj_metacell[['metacell_purity']]<0.9] = colnames(obj_metacell)[obj_metacell[['metacell_purity']]<0.9]
  }
  
  if (dub_mc_test.label) {
    if ('mcRigor' %in% colnames(obj_metacell@meta.data)) {
      centroids[['dub_mc_test']] = NA
      centroids[['dub_mc_test']][obj_metacell$mcRigor == 'dubious'] = colnames(obj_metacell)[obj_metacell$mcRigor == 'dubious']
    } else stop('Cannot label dubious metacells since test results are not provided.')
  }
  
  if(is.null(color_field)) {
    color_field <- "all"
    centroids[[color_field]] <- "lightgrey"
    
    scCoord <- data.frame(scCoord)
    scCoord[[color_field]] <- 'lightgrey'
  } else {
    if (!(color_field %in% names(obj_metacell[[]]))) stop('The field for coloring does not exist.')
    centroids[[color_field]] <- obj_metacell[[color_field]][,1]
    
    scCoord <- data.frame(scCoord)
    scCoord[[color_field]] <- obj_singlecell[[color_field]][,1]
    
    scCoord[[color_field]] = factor(scCoord[[color_field]])
    centroids[[color_field]] = factor(centroids[[color_field]], levels = levels(scCoord[[color_field]]))
  }
  
  
  
  p <- ggplot2::ggplot(scCoord,
                       ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                           colnames(scCoord)[dims[2]],
                                           color = color_field)) +
    ggplot2::geom_point(size=pt_size, alpha = sc.alpha) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))
  
  plegend = ggpubr::get_legend(ggplot2::ggplot(scCoord,
                                               ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                                                   colnames(scCoord)[dims[2]],
                                                                   color = color_field)) +
                                 ggplot2::geom_point(size=pt_size, alpha = 1) +
                                 ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2))))
  
  if (!is.null(color_field) && !is.null(cpalette)) {
    plegend = ggpubr::get_legend(ggplot2::ggplot(scCoord,
                                                 ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                                                     colnames(scCoord)[dims[2]],
                                                                     color = color_field)) +
                                   ggplot2::geom_point(size=pt_size, alpha = 1) +
                                   ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2))) + 
                                   ggplot2::scale_color_manual(values = cpalette))
  }
  
  
  
  if(continuous_metric){
    p <- p + 
      ggplot2::geom_point(data=centroids,
                          ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                              colnames(centroids)[1 + dims[2]],
                                              fill = color_field, size = metric),
                          pch=21, color = 'black', stroke = 1, alpha = mc.alpha) +
      ggplot2::scale_size_continuous(limits = c(1,max_mcsize), range = c(1,12))
    
    # p <- p + 
    #   ggplot2::geom_point(data=centroids,
    #                              ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
    #                                                  colnames(centroids)[1 + dims[2]],
    #                                                  fill = color_field, size = metric),
    #                              pch=16, alpha = mc.alpha) +
    #   ggplot2::geom_point(data=centroids,
    #                         ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
    #                                             colnames(centroids)[1 + dims[2]],
    #                                             size = metric), 
    #                         pch=1, color = 'black', stroke = 2) + 
    #   ggplot2::scale_size_continuous(limits = c(1,max_mcsize), range = c(1,12))
  }else{
    p <- p + ggplot2::geom_point(data=centroids, 
                                 ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                     colnames(centroids)[1 + dims[2]],
                                                     fill = metric), 
                                 pch=21, color = 'black', stroke = 1, alpha = mc.alpha) 
  } 
  
  if (dub_mc.label){
    p <- p + 
      ggplot2::geom_point(data=centroids[!is.na(centroids$dub_mc),],
                          ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                              colnames(centroids)[1 + dims[2]],
                                              fill = color_field, size = metric),
                          pch=21, color = 'red', stroke = 1.5, alpha = 1) +
      ggplot2::scale_size_continuous(limits = c(1,max_mcsize), range = c(1,12))
    # p <- p + ggplot2::geom_text(data = centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
    #                                                                   colnames(centroids)[1 + dims[2]],
    #                                                                   label = 'dub_mc'), color = 'black', size = 2)
    
    if(label_text){
      p <- p + ggplot2::geom_text(data = centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                                        colnames(centroids)[1 + dims[2]],
                                                                        label = 'dub_mc'), color = 'black', size = 2)
    }
  }
  
  if (dub_mc_test.label){
    p <- p + 
      ggplot2::geom_point(data=centroids[!is.na(centroids$dub_mc_test),],
                          ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                              colnames(centroids)[1 + dims[2]],
                                              fill = color_field, size = metric),
                          pch=21, color = 'red', stroke = 1.5, alpha = 1) +
      ggplot2::scale_size_continuous(limits = c(1,max_mcsize), range = c(1,12))
    
    if(label_text){
      p <- p + ggplot2::geom_text(data = centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                                        colnames(centroids)[1 + dims[2]],
                                                                        label = 'dub_mc_test'), color = 'black', size = 2)
    }
  }
  
  
  if (!is.null(color_field) && !is.null(cpalette)) {
    p <- p + 
      ggplot2::scale_fill_manual(values = cpalette) +
      ggplot2::scale_color_manual(values = cpalette) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'none', panel.grid = ggplot2::element_blank())
  } else {
    p <- p +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'none', panel.grid = ggplot2::element_blank())
  }
  
  if (!axis_lab) {
    p <- p + ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank())
  }
  
  return(list(plot = p, legend = plegend))
}




