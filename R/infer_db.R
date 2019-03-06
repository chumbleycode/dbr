#' infer_db
#'
#' The main workhorse for differential binding (db) analysis.
#'
#' @param ttT is a tidy "topTable" object, output from tidy_topTable
#' @param ttT_sub a subset - of rows - of the above, considered interesting in
#'   some way.
#' @param which_matrix the TFBM matrix: "utr1", "exonic1_utr1", "exonic1"
#' @param which_tfbms a character vector of interesting binding motifs (cols of
#'   which matrix)
#' @param n_sim the number of samples for the monte carlo inference
#' @param explicit_zeros If TRUE, then create explicit zero counts in the tfbm
#'   matrix for all genes in the sample frame that are not in tfbm database. If
#'   FALSE, then only genes the set of genes with at least one binding motif for
#'   at least one regulator.
#'
#' @return blah
#' @export
#'
#' @examples
infer_db =
  function(ttT, ttT_sub = NULL,
           which_matrix = NULL,
           which_tfbms = NULL,
           explicit_zeros = FALSE,
           n_sim = 100000){

    if(is.null(which_matrix)) which_matrix = utr1 # default matrix, if unspecified
    R   = get_matrix(which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros)
    ttT = append_matrix(ttT = ttT, ttT_sub = ttT_sub, R = R)

    # TELIS
    telis = NULL
    if(!is.null(ttT_sub)) telis <- get_telis(R = R, ttT_sub = ttT_sub, n_sim = n_sim) # possibly update telis = NULL

    # REGRESSIONS
    m =
      ttT %>%
      select(B,t,AveExpr,logFC) %>%     # THE VARIOUS OUTCOMES y
      map(regress_db,                  #
          X = select(ttT, colnames(R))) # X = THE RHS OF THE REGRESSIONS

    return(out = list(ttT   = ttT,
                      telis = telis,
                      m     = m))

    #######################################################
    # OTHER VERY ROUGH CODE: SPECULATIONS
    #######################################################

    if(0){
      t_val = results %>% select(gene, t)
      # ensure that rownames or R match names of t_value vector:
      t_value = t_val$t
      names(t_value) = t_val$gene
      # checks
      all.equal(names(t_value[rownames(R)]),  rownames(R)  )
      n_gen_tot == length(t_value)

      # null statistics
      sims = rerun(n_sim, sample(t_value))
      sims = matrix(unlist(sims), nrow = n_sim, byrow = T)
      S    = (sims / n_gen_tot) %*% R            # mean motif-count-per-gene statistic is a linear function of omega
      obs  = (t(t_value / n_gen_tot) %*% R)     # corresponding motif statistics

      # p values
      anther = NULL
      anther$p_vals_left_tail  = map2_infer_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y >= .x)) # down regulation?
      anther$p_vals_right_tail = map2_infer_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y <= .x)) # up regulation?
      anther$ttT_sub = ttT_sub
    }

    if(0){

      ########################################################
      # ANOTHER
      ########################################################

      # ensure that rownames or R match names of t_value vector:
      post_p_value = post_p_val$probability_of_effect
      names(post_p_value) = post_p_val$gene
      # checks
      all.equal(names(post_p_value[rownames(R)]),  rownames(R)  )
      n_gen_tot == length(post_p_value)

      # null statistics
      sims = rerun(n_sim, sample(post_p_value))
      sims = matrix(unlist(sims), nrow = n_sim, byrow = T)
      S    = (sims / n_gen_tot) %*% R            # mean motif-count-per-gene statistic is a linear function of omega
      obs  = (t(post_p_value / n_gen_tot) %*% R)     # corresponding motif statistics

      # p values
      p_vals_left_tail  = map2_infer_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y >= .x)) # down regulation?
      p_vals_right_tail = map2_infer_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y <= .x)) # up regulation?


      voom(counts = counts[, complete.cases(design)],
           design = design[complete.cases(design), ]) %>% # arrayWeights %>%
        limma::lmFit %>%
        eBayes %>%
        topTable(coef = of_in, n = Inf) %>%
        as_tibble(rownames = "gene") %>%
        gf_histogram(~exp(B)/(1+exp(B))) %>%
        gf_labs(title = "probability of an effect")

    }

    ########################################################
    # MORE INDULGENT STATISTICAL EXPLORATIONS
    ########################################################

    if(0){
      par(mfrow = c(3, 3))
      for(ii in union(p_vals_left_tail[p_vals_left_tail <= 0.05/2] %>% names, p_vals_right_tail[p_vals_right_tail <= 0.05/2] %>% names) ){
        S[, ii] %>% hist(100, main = ii)
        abline(v = obs[colnames(obs) == ii])
      }
    }


    if(0){

      # NULL DISTRIBUTION
      S1 = S[1:(n_sim/2), ]
      S2 = S[(n_sim/2+1):n_sim, ]
      null_p_vals = vector("list", 100)
      null_dn     = vector("list", 100)
      for(ii in 1:100){
        null_obs = S2[ii, ]
        null_p_vals[[ii]] = map2_infer_dbl(.x = as_tibble(S1), .y = null_obs, ~ mean(.y <= .x))
        null_dn[[ii]] = S1[, ii] %>% unlist
      }
      (null_dn[[1]]*n) %>% hist(1000, main = "total tfbm ALX3 count in df gene set")
      null_dn[[1]] %>% hist(1000,  main = "average tfbm ALX3 count per df gene")
      null_p_vals %>% map(1) %>% unlist %>% hist(1000,  main = "average tfbm ALX3 count per df gene")

      null_p_vals[1] %>% unlist %>% hist(1000) # a single experiment over all
      null_p_vals  %>% unlist %>% hist(main = "p-value null distribution over")

      # Steve's analytic approximation to the permutation distribution
      diff_between_analytic_and_permutation_mean = NULL
      diff_between_analytic_and_permutation_var = NULL
      for(ii in 1:dim(S)[2]){
        rbind(c(mean(S[, ii]), sd(S[, ii])) ,
              c(mean(R[, ii]), sd(R[, ii])/sqrt(n)))

        diff_between_analytic_and_permutation_mean[ii] = ( mean(S[, ii]) - mean(R[, ii]))
        diff_between_analytic_and_permutation_var[ii] = (sd(S[, ii])    - sd(R[, ii])/sqrt(n))
      }
      hist(diff_between_analytic_and_permutation_mean, xlim = c(-1,1))
      hist(diff_between_analytic_and_permutation_var, xlim = c(-1,1))

    }

    if(0){

      ind = p_vals < 0.05
      ind[40:length(ind)] = FALSE

      # mean (same same)
      S[, ind] %>%
        as_tibble %>%
        gather(tf, freq_total_tfbm_in_n_differentially_expressed_genes) %>%
        left_join(obs_tib, by = "tf") %>%
        gf_histogram(~freq_total_tfbm_in_n_differentially_expressed_genes, bins = 100) %>%
        gf_vline(xintercept = ~V1/n) %>%
        gf_facet_wrap(~tf)

    }
  }


#' append_db
#'
#' @param ttT see infer_db
#' @param ttT_sub  see infer_db
#' @param which_matrix see infer_db
#' @param which_tfbms see infer_db
#' @param explicit_zeros see infer_db
#'
#' @return
#' @export
#'
#' @examples
append_db =
  function(ttT, ttT_sub = NULL,
           which_matrix = NULL,
           which_tfbms = NULL,
           explicit_zeros = FALSE){

    if(is.null(which_matrix)) which_matrix = utr1 # default matrix, if unspecified
    R   = get_matrix(which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros)
    ttT = append_matrix(ttT = ttT, ttT_sub = ttT_sub, R = R)
  }







########################################################
# UTILITY TO BE CALLED ONLY FROM ABOVE (else resolve ...)
########################################################

#' Unexported function to LOAD AND TIDY TFBM MATRIX
#'
#' @param which_matrix blah
#' @param which_tfbms  blah
#' @param explicit_zeros blah
#'
#' @return a tidied tfbm matrix
#'
#' @examples
get_matrix <- function(which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros){

  # LOAD TFBM MATRIX
  R = which_matrix # the gene x motif binding "R"egulation matrix

  # ATTEMPT TO REMOVE UNINFORMATIVE TFBMS. REMOVE ROWS WITH NO VARIATION.
  # COULD ALSO ADDRESS COLINEARITY HERE.
  R = R[, matrixStats::colSds(R) >= 0.1]

  # PERHAPS SUBSET COLS OF R
  # CASE 1: ONLY ONE NON NULL TFBM SPECIFIED BEWARE TYPE CHANGE
  # CASE 2: ANY OTHER NON-NULL, NON-SINGLETON PROPER SUBSET OF COLS OF R
  # CASE 3: EXAMINE ALL TFBMS
  if((!is.null(which_tfbms)) & (length(which_tfbms) == 1)){
    Rt  = R[rownames(R) %in% ttT$gene, which_tfbms]
    R   = matrix(Rt, ncol = 1) %>% `rownames<-`(names(Rt)) %>% `colnames<-`(which_tfbms)
  } else if((!is.null(which_tfbms)) & (length(which_tfbms) >= 2)) {
    R = R[rownames(R) %in% ttT$gene, which_tfbms]        # restrict R to those genes which could in principle be differentially expressed
  } else if(is.null(which_tfbms)) {
    which_tfbms = colnames(R)
    R = R[rownames(R) %in% ttT$gene, which_tfbms]
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

#' Join TFBM data to ttT
#'
#' @param ttT blah
#' @param ttT_sub blah
#' @param R blah
#'
#' @return a new gene wise tidy table with appended tfbm counts
#'
append_matrix <- function(ttT = ttT, ttT_sub = ttT_sub, R = R){
  ttT =
    as_tibble(R, rownames = "gene") %>%
    left_join(ttT, by = "gene") %>%
    mutate(a_priori = gene %in% ttT_sub$gene) # WHICH WERE "INTERESTING", POSSIBLY NULL
}


########################################################
# UTILITY TO BE CALLED ONLY FROM ABOVE (else resolve ...)
########################################################


#' Title
#'
#' @param R blah
#' @param ttT_sub bla
#' @param n_sim blah
#'
#' @return blah
#'
#' @examples
get_telis <- function(R = R, ttT_sub = ttT_sub, n_sim = n_sim){

  telis = NULL
  ########################################################
  # TELIS p values
  ########################################################

  ########################################################
  # A NEW NON-PARAMETRIC SCHEME:
  # NON PARAMETRIC MONTE CARLO NULL DISTRIBUTION
  ########################################################

  responsive_lgl = rownames(R) %in% ttT_sub$gene # df gene set
  (n_gene = sum(responsive_lgl))                 # size of focal gene set within R

  sims = rerun(n_sim, sample(responsive_lgl))
  sims = matrix(unlist(sims), nrow = n_sim, byrow = T)
  sims = unique(sims) # for purists: this ensures that we only sample from the set of all possible permutations, WITH replacement

  S    = (sims / n_gene) %*% R              # mean motif-count-per-gene statistic is a linear function of omega
  obs  = (t(responsive_lgl / n_gene) %*% R) # corresponding motif statistics

  telis$npar$p_vals_left_tail  = map2_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y >= .x)) # down regulation?
  telis$npar$p_vals_right_tail = map2_dbl(.x = as_tibble(S), .y = obs, ~ mean(.y <= .x)) # up regulation?

  ########################################################
  # STEVE'S ORIGINAL PARAMETRIC IID APPROXIMATION
  ########################################################

  mu = colMeans(R)
  se = matrixStats::colSds(R)/sqrt(n_gene)
  telis$par$p_vals_left_tail  = pnorm(obs, mu, se, lower.tail = T) %>% as.vector %>% `names<-`(names(mu)) # downregulation
  telis$par$p_vals_right_tail = pnorm(obs, mu, se, lower.tail = F) %>% as.vector %>% `names<-`(names(mu)) # upregulation

  ########################################################
  # SOME REPORTING
  ########################################################

  print(str_c("dimensions of TFBM matrix R  :  ", dim(R)[1], " x ", dim(R)[2]))
  print(str_c("n genes within R             : ",  length(responsive_lgl)))
  print(str_c("n interesting genes within R : ",  sum(responsive_lgl)))

  return(telis = telis)
}
