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
#' @importFrom purrr map map_df
#' @importFrom dplyr filter mutate
#'
#' @examples
regress_db =
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
      stats::lm(y ~ ., data = X) %>%
      broom::tidy() %>%
      filter(term != "(Intercept)") # all together

    return(out = out)
  }
