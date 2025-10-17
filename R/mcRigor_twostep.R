

#' @title Step 1 of mcRigor two-step: Identifying the single cells that constitute dubious metacells and need to be re-partitioned.
#'
#' @description Step 1 of mcRigor two-step: Identifying the single cells that constitute dubious metacells and need to be re-partitioned.
#'
#' @param obj_singlecell A Seurat object of the original single cells.
#' @param sc_membership A named vector (or dataframe) of the metacell membership of single cells.
#' @param TabMC A dataframe containing the permutation results and mcRigor results in dubious metacell detection. Saved in the previous steps through function mcRigor_DETECT.
#' @param method The metacell partitioning method used to build metacells. Available: "seacells" (default, SEACells), "mc1" (MetaCell), "mc2" (MetaCell2), "supercell" (SuperCell), "metaq" (MetaQ).
#' 
#' @return A list containing the following fields: 
#' 
#' \item{obj_sc_dub}{A Seurat object of single cells that are marked as dubious and will be re-partitioned.}
#' \item{target_res}{A list containing the permutation results, derived thresholds, and mcRigor results in dubious metacell detection, given by function mcRigor_threshold.}
#' \item{obj_metacell_step1}{The metacell object built under the membership sc_membership, with more dubious metacells marked.}
#'
#' @export


mcRigorTS_Step1 <- function(obj_singlecell,
                            sc_membership,
                            TabMC = NULL,
                            method = c('seacells', 'mc1', 'mc2', 'supercell', 'metaq')) {
  
  method = match.arg(method)
  
  if (is.null(TabMC)) {
    cell_membership_temp = data.frame(sc_membership,row.names = names(sc_membership))
    names(cell_membership_temp)[1] = sub("mc(\\d+)-.*", "\\1", sc_membership[1])
    TabMC = mcRigor_DETECT(obj_singlecell = obj_singlecell, cell_membership = cell_membership_temp)$TabMC
  }
  
  mmc_res_org = mcRigor_threshold(TabMC)
  Threshold = mmc_res_org$threshold
  
  obj_metacell_step1 = mcRigor_buildmc(obj_singlecell = obj_singlecell, sc_membership = sc_membership,
                                       add_testres = T, test_stats = TabMC, test_cutoff = 0.15)
  # table(obj_metacell_bmcite@misc$cell_membership$testres_sc)
  dub_sc = rownames(obj_metacell_step1@misc$cell_membership)[obj_metacell_step1@misc$cell_membership$mcRigor_sc == 'dubious']
  
  obj_sc_dub = obj_singlecell[, match(dub_sc, colnames(obj_singlecell))]
  
  obj_sc_dub = obj_sc_dub[rowSums(GetAssayData(obj_sc_dub, layer = 'counts')) > 0, ]
  
  saveRDS(obj_sc_dub, file = paste0('dubsc_', method, '.rds'))
  
  if (method %in% c('seacells', 'mc2', 'metaq')) {
    
    if (!dir.exists(method)) dir.create(method)
    
    write.csv(t(GetAssayData(obj_sc_dub, layer = 'counts')), file = paste0(method, '/counts.csv'))
    write.csv(obj_sc_dub@meta.data, file = paste0(method, '/metadata.csv'))
    
  }
  
  return(list(obj_sc_dub = obj_sc_dub,
              target_res = mcRigor_threshold(TabMC),
              obj_metacell_step1 = obj_metacell_step1))
}





#' @title Step 2 of mcRigor two-step: Re-partition single cells that previously formed dubious metacells into smaller metacells.
#'
#' @description Step 2 of mcRigor two-step: Re-partition single cells that previously formed dubious metacells into smaller metacells.
#'
#' @param step1_res A list of results obtained form Step 1 of mcRigor two-step, given by function mcRigorTS_Step1.
#' @param obj_singlecell A Seurat object of the original single cells.
#' @param cell_membership_twostep A dataframe, each column of which contains the metacell membership of single cells that need to be re-partitioned under a given gamma (granularity level). 
#' The column names should be the corresponding gamma's (in the character type). The row names should be the single cell id's.
#' @param TabMC A dataframe containing the permutation results and mcRigor results in dubious metacell detection. Saved in the previous steps through function mcRigor_DETECT.
#' @param method The metacell partitioning method used to build metacells. Available: "seacells" (default, SEACells), "mc1" (MetaCell), "mc2" (MetaCell2), "supercell" (SuperCell), "metaq" (MetaQ).
#' @param custom_thre A boolean indicating whether to re-compute the threshold values.
#' @param color_field The variable based on which to color the cells.
#' 
#' @return A list containing the following fields: 
#' 
#' \item{obj_sc_dub}{A Seurat object of single cells that are marked as dubious and will be re-partitioned.}
#' \item{target_res}{A list containing the permutation results, derived thresholds, and mcRigor results in dubious metacell detection, given by function mcRigor_threshold.}
#' \item{obj_metacell_step1}{The metacell object built under the membership sc_membership, with more dubious metacells marked.}
#'
#' @export


mcRigorTS_Step2 <- function(step1_res, 
                            obj_singlecell,
                            cell_membership_twostep,
                            twostep_gamma = NULL,
                            method = c('seacells', 'mc1', 'mc2', 'supercell', 'metaq'), 
                            custom_thre = F,
                            color_field = NULL){
  
  method = match.arg(method)
  
  temp_res = mcRigor_DETECT(obj_singlecell = step1_res$obj_sc_dub,
                            cell_membership = cell_membership_twostep,
                            output_file = paste0(method, '_tabmc_twostep.rds'),
                            assay_type = 'RNA')
  TabMC_twostep = temp_res$TabMC
  
  if (!custom_thre) {
    mmc_res_twostep = mcRigor_tradeoff(TabMC_twostep, threshold = step1_res$target_res$threshold)
  } else {
    mmc_res_twostep = mcRigor_tradeoff(TabMC_twostep, threshold = NULL)
  }
  
  if (is.null(twostep_gamma)) twostep_gamma = mmc_res_twostep$optimized$gamma
  sc_membership_twostep = cell_membership_twostep[[as.character(twostep_gamma)]]
  names(sc_membership_twostep) = rownames(cell_membership_twostep)
  
  sc_membership_tab = step1_res$obj_metacell_step1@misc$cell_membership
  sc_membership_first = sc_membership_tab$Metacell
  names(sc_membership_first) = rownames(sc_membership_tab)
  sc_membership_first = sc_membership_first[sc_membership_tab$mcRigor_sc == 'trustworthy']
  
  TabMC = step1_res$target_res$TabMC
  tgamma = sub("mc(\\d+)-.*", "\\1", Cells(step1_res$obj_metacell_step1)[1])
  TabMC = TabMC[TabMC$gamma == as.integer(tgamma),]
  
  
  if (!custom_thre){
    
    col_names_tab = intersect(colnames(TabMC), colnames(TabMC_twostep))
    TabMC = TabMC[, col_names_tab]
    TabMC_twostep = TabMC_twostep[, col_names_tab]
    
    sc_membership_final = c(sc_membership_first, sc_membership_twostep)
    
    obj_metacell_final = mcRigor_buildmc(obj_singlecell = obj_singlecell,
                                         sc_membership = sc_membership_final,
                                         add_testres = T, test_stats = rbind(TabMC, TabMC_twostep), Thre = step1_res$target_res$threshold)
    
    plot = mcRigor_projection(obj_singlecell = obj_singlecell,
                              sc_membership = sc_membership_final,
                              add_testres = T, test_stats = rbind(TabMC, TabMC_twostep), Thre = step1_res$target_res$threshold,
                              dub_mc_test.label = T,
                              color_field = color_field)
    
  } else {
    
    testres_1 = rep('trustworthy', length(unique(sc_membership_first)))
    names(testres_1) = unique(sc_membership_first)
    testres_2 = mmc_res_twostep$TabMC$mcRigor
    names(testres_2) = rownames(mmc_res_twostep$TabMC)
    testres_2 = testres_2[names(testres_2) %in% unique(sc_membership_twostep)]
    
    testres_final = c(testres_1, testres_2)
    
    sc_membership_final = c(sc_membership_first, sc_membership_twostep)
    
    obj_metacell_final =  mcRigor_buildmc(obj_singlecell = obj_singlecell,
                                          sc_membership = sc_membership_final,
                                          add_testres = T, testres = testres_final)
    
    plot = mcRigor_projection(obj_singlecell = obj_singlecell,
                              sc_membership = sc_membership_final,
                              add_testres = T, testres = testres_final,
                              dub_mc_test.label = T,
                              color_field = color_field)
    
  }
  
  mc_sc_id = rownames(obj_metacell_final@misc$cell_membership)[obj_metacell_final@misc$cell_membership$mcRigor_sc == 'dubious']
  sc_membership_twostep_withsc = sc_membership_twostep
  sc_membership_twostep_withsc[mc_sc_id] = mc_sc_id
  
  sc_membership_final_withsc = c(sc_membership_first, sc_membership_twostep_withsc)
  
  obj_metacell_final_withsc = mcRigor_buildmc(obj_singlecell = obj_singlecell,
                                              sc_membership = sc_membership_final_withsc,
                                              add_testres = F)
  
  plot_withsc = mcRigor_projection(obj_singlecell = obj_singlecell,
                                   sc_membership = sc_membership_final_withsc,
                                   add_testres = F,  dub_mc_test.label = F,
                                   color_field = color_field)
  
  
  return(list(obj_metacell_final = obj_metacell_final,
              obj_metacell_final_withsc = obj_metacell_final_withsc,
              gamma1 = tgamma,
              gamma2 = twostep_gamma,
              mmc_res_twostep = mmc_res_twostep,
              plot = plot,
              plot_sc = plot_withsc))
}
