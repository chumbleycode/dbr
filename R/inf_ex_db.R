
#' inf_ex_db
#'
#' A convenience function
#'
#' Used to create a pipeline around infer_db.
#'
#' @param rhs the right hand side of a regression, as a character vector of variable names
#' @param of_in  see infer_db
#' @param which_matrix  see infer_db
#' @param dat  an expression set object
#' @param n_sim  see infer_db
#' @param ... anything else
#'
#' @return blah
#' @export
#'
#' @import broom
#' @importFrom dplyr %>%
#' @importFrom Biobase exprs pData
#' @importFrom limma voom eBayes
#' @importFrom stringr str_c str_detect
#' @examples
inf_ex_db =
  function(rhs, of_in = of_in, which_matrix = utr1, dat, n_sim = 10000, ...){

    # FILTER OUT INADMISSIBLE GENES, DEFINE OUTCOME AND DESIGN
    e_genes = edgeR::filterByExpr(Biobase::exprs(dat))
    e_genes = names(e_genes[e_genes == TRUE])
    counts  = dat[e_genes, ] %>% Biobase::exprs()
    phen    = dat[e_genes, ] %>% Biobase::pData()
    design  = stats::model.matrix(stats::as.formula(stringr::str_c("~", rhs)), data = phen)

    # ROBUSTNESS TO MINOR DEVIATIONS IN PARAMETER NAME PRODUCED BY LMFIT BELOW, E.G. PSMOKING -> PSMOKING1
    of_in   =  colnames(design)[which(stringr::str_detect(colnames(design),unlist(of_in)))]

    # ERROR HANDLING: IF DESIGN[, OF_IN] IS A MATRIX, THEN OF_IN INDEXES MORE THAN ONE COLLUMN OF DESIGN.
    if (is.matrix(design[, of_in])) stop("variable of interest (of_in) refers to multiple columns of the design!")
    print(of_in)

    # GENERAL RESULTS
    ok = stats::complete.cases(design)
    ttT =
      voom(counts = counts[, ok],
           design = design[ok, ]) %>% # arrayWeights %>%
      limma::lmFit() %>%
      limma::eBayes() %>%
      limma::topTable(coef = of_in, number = Inf) %>%
      dplyr::as_tibble(rownames = "gene")

    # INFER DIFFERENTIAL TFBM
    out =
      ttT %>%
      infer_db(ttT_sub      = dplyr::filter(ttT, P.Value <= 0.05),
               which_matrix = which_matrix,
               n_sim        = 10000,
               ...) %>%
      extract_db()

    return(out = out)

  }
