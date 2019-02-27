
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
