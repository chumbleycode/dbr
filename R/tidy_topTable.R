
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
  x %>%
    limma::topTable(coef = of_in, n = Inf, ...) %>%
    dplyr::as_tibble(rownames = "gene")
}
