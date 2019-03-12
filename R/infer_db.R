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
#' @param perm_telis if ttT_sub is provided, do you want non-parameteric permutation TeLiS in addition to parametric
#'
#' @return blah
#'
#' @export
#'
#' @examples
infer_db =
  function(ttT, ttT_sub = NULL,
           which_matrix = NULL,
           which_tfbms = NULL,
           explicit_zeros = TRUE,
           perm_telis = FALSE,
           n_sim = 100000){

    if(is.null(which_matrix)) which_matrix = utr1 # default matrix, if unspecified
    R   = get_matrix(ttT = ttT, which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros)
    ttT = append_matrix(ttT = ttT, ttT_sub = ttT_sub, R = R)

    # TELIS
    telis = NULL
    if(!is.null(ttT_sub)) telis <- get_telis(R = R, ttT_sub = ttT_sub, n_sim = n_sim, perm_telis = perm_telis) # possibly update telis = NULL

    # REGRESSIONS
    m =
      ttT %>%
      dplyr::select(B,t,AveExpr,logFC) %>%     # THE VARIOUS OUTCOMES y
      purrr::map(regress_db,                  #
          X = dplyr::select(ttT, colnames(R))) # X = THE RHS OF THE REGRESSIONS

    return(out = list(ttT   = ttT,
                      telis = telis,
                      m     = m))
  }


