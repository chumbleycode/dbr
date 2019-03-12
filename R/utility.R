
########################################################
# PRETTY WAY TO TABULATE TFBM PVALUES
########################################################

#' Title
#'
#' @param pval   blah
#'
#' @return  blah
#'
#' @examples
disp =
  function(pval){
    # simulation and multiple correction parameters
    method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",   "fdr", "none")
    method = method %>% `names<-`(method)
    method %>%
      purrr::map(~p.adjust(pval, method = .x, n = length(pval))) %>%
      unlist %>%
      `[`(. <= 0.05) %>%
      sort %>%
      tibble::enframe("tf", "p_value") %>%
      knitr::kable()
  }





########################################################
# UTILITY TO BE CALLED ONLY FROM ABOVE (else resolve ...)
########################################################

#' append_matrix
#'
#' Join TFBM data to ttT
#'
#' @param ttT a tidied topTable object
#' @param ttT_sub a tidied topTable object whose rows contain only "DE" genes
#' @param R a gene by binding motif count matrix
#'
#' @return a new gene wise tidy table with appended tfbm counts
#'
append_matrix <- function(ttT = ttT, ttT_sub = ttT_sub, R = R){

  ttT =
    as_tibble(R, rownames = "gene") %>%
    left_join(ttT, by = "gene") %>%
    mutate(gene_subset = gene %in% ttT_sub$gene) %>% # A la TeLiS, possibly null
    dplyr::select(names(ttT), gene_subset, colnames(R)) # reorder collumns

  return(ttT = ttT)
}

########################################################
# UTILITY TO BE CALLED ONLY FROM ABOVE (else resolve ...)
########################################################

#' get_matrix
#'
#' Unexported function to LOAD AND TIDY TFBM MATRIX
#'
#' @param which_matrix which gene by binding motif count matrix
#' @param which_tfbms  named columns of "which_matrix"
#' @param explicit_zeros Add explicit zero counts to the tfbm matrix for genes in sampling frame, but not in tfbm matrix
#' @param ttT a tidied topTable object
#'
#' @return a tidied tfbm matrix
#'
#' @examples
get_matrix <- function(ttT = ttT, which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros){

  # LOAD TFBM MATRIX
  R = which_matrix # the gene x motif binding "R"egulation matrix

  # ATTEMPT TO REMOVE UNINFORMATIVE TFBMS. REMOVE ROWS WITH NO VARIATION.
  # COULD ALSO ADDRESS COLINEARITY HERE.
  R = R[, matrixStats::colSds(R) >= 0.1]

  # PERHAPS SUBSET COLS OF R
  # CASE 1: ONLY ONE NON NULL TFBM SPECIFIED BEWARE TYPE CHANGE
  # CASE 2: ANY OTHER NON-NULL, NON-SINGLETON PROPER SUBSET OF COLS OF R
  # CASE 3: EXAMINE ALL TFBMS
  in_frame  = rownames(R) %in% ttT$gene # which tfbm genes are in the sampling frame
  if((!is.null(which_tfbms)) & (length(which_tfbms) == 1)){
    Rt  = R[in_frame, which_tfbms]
    R   = matrix(Rt, ncol = 1) %>% `rownames<-`(names(Rt)) %>% `colnames<-`(which_tfbms)
  } else if((!is.null(which_tfbms)) & (length(which_tfbms) >= 2)) {
    R = R[in_frame, which_tfbms]        # restrict R to those genes which could in principle be differentially expressed
  } else if(is.null(which_tfbms)) {
    which_tfbms = colnames(R)
    R = R[in_frame, which_tfbms]
  }

  # The original base tfbm matrix only includes genes with at least one tfbm for at least one regulator.
  # Genes not mentioned in the base tfbm matrix implicitly have a count of zero for every tfbm.
  # The preceding code additionally excludes any genes in the base tfbm matrix but not in the sampling frame (ttT$gene).
  # The next option introduces explicit zeros to R, for all genes in the sampling frame but with (implicitly) zero tfbm count for any motif.
  if(explicit_zeros){
    print("Adding explicit zero counts to the tfbm matrix for genes in sampling frame, but not in tfbm matrix")
    not_in_R = setdiff(ttT$gene, rownames(R))
    nR = matrix(0,length(not_in_R), dim(R)[2])
    rownames(nR) = not_in_R
    R = rbind(R, nR)
  }
  return(R = R)
}


########################################################
# UTILITY TO BE CALLED ONLY FROM ABOVE (else resolve ...)
########################################################

#' get_telis
#'
#' @param R a gene by binding motif count matrix
#' @param ttT_sub a tidied topTable object whose rows contain only "DE" genes
#' @param n_sim number of permutations if non parametric telis
#' @param perm_telis non-parametric telis?
#'
#'
#' @examples
get_telis <- function(R = R, ttT_sub = ttT_sub, n_sim = n_sim, perm_telis = FALSE){

  telis = NULL

  ########################################################
  # TELIS p values
  ########################################################

  responsive_lgl = rownames(R) %in% ttT_sub$gene # df gene set
  (n_gene = sum(responsive_lgl))                 # size of focal gene set within R
  obs  = (t(responsive_lgl / n_gene) %*% R) # corresponding motif statistics

  ########################################################
  # STEVE'S ORIGINAL PARAMETRIC IID TELIS (THE GAUSSIAN )
  ########################################################

  mu = colMeans(R)
  se = matrixStats::colSds(R)/sqrt(n_gene)
  telis$par$p_under = stats::pnorm(obs, mu, se, lower.tail = T) %>% as.vector %>% `names<-`(names(mu)) # downregulation
  telis$par$p_over  = stats::pnorm(obs, mu, se, lower.tail = F) %>% as.vector %>% `names<-`(names(mu)) # upregulation

  if(perm_telis) {

    if (dim(R)[1] > 100) {
      print("For large sample frames, there is little benefit in return for the computational expense of permutation analysis")
    }
    ########################################################
    # A NEW NON-PARAMETRIC SCHEME:
    # NON PARAMETRIC MONTE CARLO NULL DISTRIBUTION
    ########################################################

    sims = purrr::rerun(n_sim, sample(responsive_lgl))
    sims = matrix(unlist(sims), nrow = n_sim, byrow = T)
    sims = unique(sims) # for purists: this ensures that we only sample from the set of all possible permutations, WITH replacement

    S    = (sims / n_gene) %*% R              # mean motif-count-per-gene statistic is a linear function of omega

    telis$npar$p_under  = purrr::map2_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y >= .x)) # down regulation?
    telis$npar$p_over   = purrr::map2_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y <= .x)) # up regulation?

  }

  ########################################################
  # SOME REPORTING
  ########################################################

  print(str_c("dimensions of TFBM matrix R  : ",  dim(R)[1], " x ", dim(R)[2]))
  print(str_c("n genes within R             : ",  length(responsive_lgl)))
  print(str_c("n interesting genes within R : ",  sum(responsive_lgl)))

  return(telis = telis)
}

