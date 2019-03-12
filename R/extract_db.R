
########################################################
# EXTRACT THE OUTPUT OF FUNCTION infer_db()
########################################################

#' extract_db
#'
#' A TABLE OF P-VALUES FOR DIFFERENT METHODS:
#'

#'
#' @param x   the output of \code{\link{infer_db}}
#' @param methods  depricated
#' @param which_outcome in general the outcome variable may be chosen as logFC, B, etc
#'
#' @return
#'
#'  A tibble whose columns give tfbm p-values for various methods:
#'
#' \itemize{
#'  \item p_par two-tailed p-value from parametric telis, returned  only if argument 'tt_sub' is given to \code{\link[dbr]{infer_db}}.
#'  \item par_p_over  = upper tailed one-sided p-value for p_par
#'  \item par_p_under = lower tailed one-sided p-value for p_par
#'  \item p_uni  = simple univariate regression of b on each tfbm
#'  \item p_cov  = multple regression of b on all tfbm's simultaniously
#'  \item p_npar = nonparametric telis, returned  only if infer_db() argument 'tt_sub' is provided and perm_telis = f
#' }
#'
#' @export
#'
#' @import broom
#' @importFrom dplyr left_join
#' @importFrom purrr is_empty
#' @importFrom tibble enframe
#'
#' @seealso \code{\link[dbr]{infer_db}}
#' @examples
extract_db =
  function(x, which_outcome = "logFC", methods = NULL){
    outcome = x$m[[which_outcome]]

    # Works even if all parameters in the multiple regression are identified (na in output, so the dimension of p_uni and p_cov are mismatched)
    out =
      outcome$m_uni %>%
      dplyr::select(tfbm = term,
             estimate_uni = estimate,
             p_uni = p.value) %>%
      dplyr::left_join(select(outcome$m_cov,
                       tfbm = term,
                       estimate_cov = estimate,
                       p_cov = p.value), by = "tfbm")

    if(!is.null(x$telis)){

      # Estimate/append telis
      out_telis =
        x$telis %>%
        unlist %>%
        tibble::enframe() %>%
        tidyr::separate(name, c("method", "side", "tfbm"), sep = "\\.") %>%
        tidyr::unite("method", c("method", "side")) %>%
        tidyr::spread(method, value) %>%
        dplyr::mutate(p_par  =  2 * pmin(par_p_under, par_p_over)) #,
               # p_npar =  NA) # check: 2 * for two tails

      if(!is_empty(x$telis$npar)){

        # If costly permutation analysis is done, overwrite NA
         out_telis =
          out_telis %>%
          dplyr::mutate(p_npar = pmin(1, 2 * pmin(npar_p_under, npar_p_over)))

      }
      out =
        out %>%
        dplyr::left_join(out_telis, by = "tfbm")

    }

    if(!is.null(methods)) out = out %>% dplyr::select(tfbm, methods) # possibly subset
    return(out = out)
  }

