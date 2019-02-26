
########################################################
# A CONVENIENCE FOR CALLING infer_db() REPEATEDLY (SEE compare_model_code.R)
########################################################

#' get_tfbm_p_vals
#'
#' Used to create a pipeline around infer_db.
#'
#' @param rhs  blah
#' @param of_in  see infer_db
#' @param which_matrix  see infer_db
#' @param dat  an expression set object
#' @param n_sim  see infer_db
#'
#' @return blah
#' @export
#'
#' @import broom
#' @importFrom dplyr %>%
#' @importFrom Biobase exprs pData
#' @importFrom limma voom eBayes
#' @examples
get_tfbm_p_vals =
  function(rhs, of_in = of_in, which_matrix = utr1, dat, n_sim = 10000){

    e_genes = edgeR::filterByExpr(Biobase::exprs(dat))
    e_genes = names(e_genes[e_genes == TRUE])

    # FILTER OUT INADMISSIBLE GENES, DEFINE OUTCOME AND DESIGN
    counts    = dat[e_genes, ] %>% Biobase::exprs()
    phen      = dat[e_genes, ] %>% Biobase::pData()
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
      limma::lmFit() %>%
      eBayes %>%
      topTable(coef = of_in, n = Inf) %>%
      as_tibble(rownames = "gene")

    ########################################################
    # INFER DIFFERENTIAL TFBM
    ########################################################

    out =
      ttT %>%
      infer_db(ttT_sub        = filter(ttT, P.Value <= 0.05),
               which_matrix   = which_matrix,
               n_sim          = 10000) %>%
      extract_db()

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
tidy_topTable = function(x, of_in, ...) {
 # browser()
  x %>%
    limma::topTable(coef = of_in, n = Inf, ...) %>%
    dplyr::as_tibble(rownames = "gene")
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
