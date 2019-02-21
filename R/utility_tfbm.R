
########################################################
# TIDY TOP TABLE OUTPUT
########################################################

#' tidy_topTable
#'
#' This simply reformats topTable's output
#'
#' @param x the output of limma::eBays
#' @param of_in # string refering to a single collumn of design matrix by name
#' @param ... any other argument to toptable
#'
#' @return a reformatted version of topTable
#' @export
#'
#' @importFrom dplyr as_tibble
#' @importFrom limma topTable
#' @examples
#'
tidy_topTable = function(x, of_in, ...) {

  x %>%
    limma::topTable(coef = of_in, n = Inf, ...) %>%
    dplyr::as_tibble(rownames = "gene")
}

########################################################
# THE MAIN WORKHORSE FOR DIFFERENTIAL BINDING (DB) ANALYSIS
########################################################

#' Title
#'
#' @param ttT is a tidy "topTable" object, output from tidy_topTable
#' @param ttT_sub a subset - of rows - of the above, considered interesting in some way.
#' @param which_matrix the TFBM matrix: "utr1", "exonic1_utr1", "exonic1"
#' @param which_tfbms a character vector of interesting binding motifs (cols of which matrix)
#' @param n_sim the number of samples for the monte carlo inference
#'
#' @return blah
#' @export
#'
#' @examples
infer_db =
  function(ttT, ttT_sub = NULL,
           which_matrix = NULL,
           which_tfbms = NULL,
           n_sim = 10000 ){

    if(is.null(which_matrix)) which_matrix = utr1 # default matrix, if unspecified

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

    # IS TELIS REQUESTED?
    telis = NULL
    if(!is.null(ttT_sub)){

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

    }

    ########################################################
    # tidy topTable
    ########################################################

    ttT =
      as_tibble(R, rownames = "gene") %>%
      left_join(ttT, by = "gene") %>%
      mutate(a_priori = gene %in% ttT_sub$gene) # WHICH WERE "INTERESTING", POSSIBLY NULL

    ########################################################
    # GENE POPULATION (OR WHOLE GENOME) SECOND LEVEL REGRESSIONS
    ########################################################

    m =
      ttT %>%
      select(B,t,AveExpr,logFC) %>%     # THE VARIOUS OUTCOMES y
      map(regress,                  #
          X = select(ttT, colnames(R))) # X = THE RHS OF THE REGRESSIONS

    return(out = list(ttT   = ttT,
                      telis = telis,
                      m     = m))


    #######################################################
    # OTHER SPECULATIONS
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



########################################################
# A CONVENIENCE FOR CALLING infer_db() REPEATEDLY (SEE compare_model_code.R)
########################################################

#' get_tfbm_p_vals
#'
#' Used to create a pipeline around infer_db
#'
#' @param rhs  blah
#' @param of_in  blah
#' @param which_matrix  blah
#' @param dat  blah
#' @param n_sim  blah
#'
#' @return blah
#' @export
#'
#' @import broom
#' @import tidyverse
#' @examples
get_tfbm_p_vals =
  function(rhs, of_in = of_in, which_matrix = "utr", dat, n_sim = 10000){

    e_genes = filterByExpr(exprs(dat))
    e_genes = names(e_genes[e_genes == TRUE])

    # FILTER OUT INADMISSIBLE GENES, DEFINE OUTCOME AND DESIGN
    counts    = dat[e_genes, ] %>% exprs
    phen      = dat[e_genes, ] %>% pData
    design    = model.matrix(as.formula(str_c("~", rhs)), data = phen)

    # ROBUSTNESS TO MINOR DEVIATIONS IN PARAMETER NAME PRODUCED BY LMFIT BELOW, E.G. PSMOKING -> PSMOKING1
    of_in    =  colnames(design)[which(str_detect(colnames(design),unlist(of_in)))]
    print(of_in)

    # ERROR HANDLING: IF DESIGN[, OF_IN] IS A MATRIX, THEN OF_IN INDEXES MORE THAN ONE COLLUMN OF DESIGN.
    if (is.matrix(design[, of_in])) {
      stop("This function cannot generally handle the case where the variable of interest (of_in) refers to multiple columns of the design!")
    }

    ########################################################
    # GENERAL RESULTS
    ########################################################

    ttT =
      voom(counts = counts[, complete.cases(design)],
           design = design[complete.cases(design), ]) %>% # arrayWeights %>%
      limma::lmFit %>%
      eBayes %>%
      topTable(coef = of_in, n = Inf) %>%
      as_tibble(rownames = "gene")

    ########################################################
    # INFER DIFFERENTIAL TFBM
    ########################################################

    out =
      ttT %>%
      infer_db(ttT_sub      =  filter(ttT, P.Value <= 0.05) ,
         which_matrix = which_matrix,
         n_sim        = n_sim)

    ########################################################
    # EXIT FUNCTION
    ########################################################

    return(out = out)

  }


########################################################
# PERFORM UNIVARIATE/MULTIVARIATE REGRESSION FOR infer_db()
########################################################

#' Title
#' Estimates simple and multiple regression of gene-specific d.e.  estimates on tfbm(s), and the regression on the total number of tfbm sites, indep of which tf they belong to.
#'
#' @param y blah
#' @param X blah
#'
#' @return blah
#' @export
#'
#' @import broom
#' @import tidyverse
#'
#' @examples
regress =
  function(y, X){
    # y = rnorm(100)
    # X = tibble(x1 = rnorm(100), x2 = rnorm(100))
    #  regress(y, X)

    out       = NULL
    out$m_uni =
      X %>%
      map(~lm(y ~ .x)) %>%
      map(tidy) %>%
      map(filter, term != "(Intercept)") %>%
      map_df(~.x) %>%
      mutate(term = names(X))
    out$m_cov =
      lm(y ~ ., data = X) %>%
      tidy %>%
      filter(term != "(Intercept)") # all together
    out$m_tot =
      lm(y ~ rowSums(X)) %>%
      tidy %>%
      filter(term != "(Intercept)") # total number of sites

    return(out = out)
  }

########################################################
# EXTRACT THE OUTPUT OF FUNCTION infer_db()
########################################################

#' Title
#'
#' @param x   blah
#' @param methods  blah
#'
#' @return  blah
#' @export
#'
#' @import broom
#' @import tidyverse
#'
#' @examples
extract_db =
  function(x, methods = NULL){
    # x is returned from infer_db()

    ########################################################
    # THE FUNCTION "extracter" GIVES A TABLE OF P-VALUES FOR DIFFERENT METHODS:
    # p_uni  = SIMPLE UNIVARIATE REGRESSION OF B ON EACH TFBM
    # p_cov  = MULTPLE REGRESSION OF B ON ALL TFBM's SIMULTANIOUSLY
    # p_tot  = SIMPLE UNIVARIATE OF B ON THE TOTAL NUMBER SITES (NOT SPECIFIC TO ONE TFBM)
    # p_par  = parametric telis, RETURNED  ONLY IF ARGUMENT 'tT_sub' IS GIVEN TO "db"
    # p_npar = nonparametric telis, RETURNED  ONLY IF ARGUMENT 'tT_sub' IS GIVEN TO "db"
    ########################################################

# NOTE: p_tot IS JUST ONE REGRESSION (OF DIFFERENTIAL EFFECT), THE RESULT IS
# REPLICATED FOR CONVENIENCE

    # print("CHECK THE 2 * below!!!!!")

    if(!is.null(x$telis)){
      out =
        x$telis %>%
        unlist %>%
        enframe %>%
        tidyr::separate(name, c("method", "side", "tfbm"), sep = "\\.") %>%
        unite("method", c("method", "side")) %>%
        spread(method, value) %>%
        mutate(p_par  = 2 * pmin(par_p_vals_left_tail, par_p_vals_right_tail),
               p_npar = 2 * pmin(npar_p_vals_left_tail, npar_p_vals_right_tail),
               p_uni  = x$m$B$m_uni$p.value,
               p_cov  = x$m$B$m_cov$p.value,
               p_tot  = x$m$B$m_tot$p.value)

    } else {

      out =
        tibble(tfbm = x$m$B$m_uni$term,
               p_uni  = x$m$B$m_uni$p.value,
               p_cov  = x$m$B$m_cov$p.value,
               p_tot  = x$m$B$m_tot$p.value)
    }

    if(!is.null(methods)) out = out %>% select(tfbm, methods) # possibly subset
    return(out = out)
  }


########################################################
# PRETTY WAY TO TABULATE TFBM PVALUES
########################################################

#' Title
#'
#' @param pval   blah
#'
#' @return  blah
#' @export
#'
#' @examples
disp =
  function(pval){
    # simulation and multiple correction parameters
    method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",   "fdr", "none")
    method = method %>% `names<-`(method)
    method %>%
      map(~p.adjust(pval, method = .x, n = length(pval))) %>%
      unlist %>%
      `[`(. <= 0.05) %>%
      sort %>%
      enframe("tf", "p_value") %>%
      knitr::kable()
  }
