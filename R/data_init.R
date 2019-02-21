
#' cleaner
#'
#' see below
#'
#' @param x a tidy object (see below)
#'
#' @return blah
#' @import tidyverse
#' @examples
cleaner = function(x){

  # A FORMATTING FUNCTION, SEE BELOW
  x =
    x %>%
    unite(tf, tf1, tf2, sep = "_", remove = F) %>%
    mutate(tf = str_replace_all(tf, "_$", "")) %>%
    dplyr::count(hgnc, tf) %>%
    spread(tf, n, fill = 0)

  x =
    x %>%
    select(-hgnc) %>%
    as.matrix %>%
    `rownames<-`(x$hgnc)
}


#' make_the_matrices
#'
#' makes all the tfbm matrices
#'
#' @return
#' @import tidyverse
#' @examples
make_the_matrices = function(){

  path_to_levitt_biomart = "../../levitt_biomart/tfbs.final.011619" # path to raw tfbm data from biomart
  tfbs0  = data.table::fread(path_to_levitt_biomart) %>% as_tibble

  ########################################################
  # THREE DIFFERENT TYPES OF MATRIX
  ########################################################

  # 1.
  utr1 =
    tfbs0 %>%
    dplyr::filter(hgnc != "", # remove dark matter
                  utr  == 1) %>%
    cleaner
  # 2.
  exonic1 =
    tfbs0 %>%
    dplyr::filter(hgnc != "", # remove dark matter
                  exonic  == 1) %>%
    cleaner
  # 3.
  exonic1_utr1 =
    tfbs0 %>%
    dplyr::filter(hgnc != "", # remove dark matter
                  exonic == 1,
                  utr  == 1) %>%
    cleaner

  ########################################################
  # SAVE
  ########################################################
  usethis::use_data(utr1, exonic1, exonic1_utr1)
  usethis::use_data(utr1, exonic1, exonic1_utr1, internal = TRUE)

  ########################################################
  # RETURN
  ########################################################

  rm(tfbs0)

}
