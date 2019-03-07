
########################################################
# EXTRACT THE OUTPUT OF FUNCTION infer_db()
########################################################

#' Title
#'
#' GIVES A TABLE OF P-VALUES FOR DIFFERENT METHODS:
#' p_uni  = SIMPLE UNIVARIATE REGRESSION OF B ON EACH TFBM
#' p_cov  = MULTPLE REGRESSION OF B ON ALL TFBM's SIMULTANIOUSLY
#' p_tot  = SIMPLE UNIVARIATE OF B ON THE TOTAL NUMBER SITES (NOT SPECIFIC TO ONE TFBM)
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
#' @import rlang
#'
#' @examples
extract_db =
  function(x, which_outcome = "B", methods = NULL){
    # x is returned from infer_db()

    # NOTE: p_tot IS JUST ONE REGRESSION (OF DIFFERENTIAL EFFECT), THE RESULT IS
    # REPLICATED FOR CONVENIENCE

    outcome = x$m[[which_outcome]]

    # THE BELOW APPEARS SIMPLER, BUT IF NOT ALL PARAMETERS IN THE MULTIPLE REGRESSION ARE IDENTIFIED (NA IN OUTPUT), THEN THE DIMENSION OF p_uni AND p_cov ARE MISMATCHED.
    # out =
    #   tibble(tfbm = outcome$m_uni$term,
    #          p_uni  = outcome$m_uni$p.value,
    #          p_cov  = outcome$m_cov$p.value,
    #          p_tot  = outcome$m_tot$p.value)
    out =
      outcome$m_uni %>%
      select(tfbm = term, estimate_uni = estimate, p_uni = p.value) %>%
      left_join(select(outcome$m_cov, tfbm = term, estimate_cov = estimate, p_cov = p.value), by = "tfbm") %>%
      mutate(p_tot  = outcome$m_tot$p.value)

    # print("CHECK THE 2 * below")
    if(!is.null(x$telis)){

      # ESTIMATE AND APPEND TELIS
      out_telis =
        x$telis %>%
        unlist %>%
        enframe %>%
        tidyr::separate(name, c("method", "side", "tfbm"), sep = "\\.") %>%
        unite("method", c("method", "side")) %>%
        spread(method, value) %>%
        mutate(p_par  =  2 * pmin(par_p_vals_left_tail, par_p_vals_right_tail),
               p_npar =  NA) # 2 * for two tails

      if(!is_empty(x$telis$npar)){

        # If costly permutation analysis is done, overwrite NA
         out_telis =
          out_telis %>%
          mutate(p_npar = pmin(1, 2 * pmin(npar_p_vals_left_tail, npar_p_vals_right_tail)))

      }
      out =
        out %>%
        left_join(out_telis, by = "tfbm")

    }

    if(!is.null(methods)) out = out %>% select(tfbm, methods) # possibly subset
    return(out = out)
  }

