
########################################################
# EXTRACT THE OUTPUT OF FUNCTION infer_db()
########################################################

#' extract_db
#'
#' A TABLE OF P-VALUES FOR DIFFERENT METHODS:
#' p_uni  = SIMPLE UNIVARIATE REGRESSION OF B ON EACH TFBM
#' p_cov  = MULTPLE REGRESSION OF B ON ALL TFBM's SIMULTANIOUSLY
#' p_par  = parametric telis, RETURNED  ONLY IF ARGUMENT 'tT_sub' IS GIVEN TO "db"
#' p_npar = nonparametric telis, RETURNED  ONLY IF ARGUMENT 'tT_sub' IS GIVEN TO "db"
#'
#' @param x   blah
#' @param methods  blah
#' @param which_outcome in general the outcome variable may be chosen as logFC, B, etc
#'
#' @return  blah
#' @export
#'
#' @import broom
#' @import tidyverse
#'
#' @examples
extract_db =
  function(x, which_outcome = "logFC", methods = NULL){
    outcome = x$m[[which_outcome]]

    # Works even if all parameters in the multiple regression are identified (na in output, so the dimension of p_uni and p_cov are mismatched)
    out =
      outcome$m_uni %>%
      select(tfbm = term,
             estimate_uni = estimate,
             p_uni = p.value) %>%
      left_join(select(outcome$m_cov,
                       tfbm = term,
                       estimate_cov = estimate,
                       p_cov = p.value), by = "tfbm")

    if(!is.null(x$telis)){

      # Estimate/append telis
      out_telis =
        x$telis %>%
        unlist %>%
        enframe %>%
        tidyr::separate(name, c("method", "side", "tfbm"), sep = "\\.") %>%
        unite("method", c("method", "side")) %>%
        spread(method, value) %>%
        mutate(p_par  =  2 * pmin(par_p_under, par_p_over)) #,
               # p_npar =  NA) # check: 2 * for two tails

      if(!is_empty(x$telis$npar)){

        # If costly permutation analysis is done, overwrite NA
         out_telis =
          out_telis %>%
          mutate(p_npar = pmin(1, 2 * pmin(npar_p_under, npar_p_over)))

      }
      out =
        out %>%
        left_join(out_telis, by = "tfbm")

    }

    if(!is.null(methods)) out = out %>% select(tfbm, methods) # possibly subset
    return(out = out)
  }

