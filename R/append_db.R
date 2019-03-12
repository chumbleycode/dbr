#' append_db
#'
#' This function takes a reformatted topTable object and appends per-gene TFBM count data.
#' @param ttT see infer_db
#' @param ttT_sub  see infer_db
#' @param which_matrix see infer_db
#' @param which_tfbms see infer_db
#' @param explicit_zeros see infer_db
#'
#' @return
#' @export
#'
#' @examples
append_db =
  function(ttT, ttT_sub = NULL,
           which_matrix = NULL,
           which_tfbms = NULL,
           explicit_zeros = TRUE){

    if(is.null(which_matrix)) which_matrix = utr1 # default matrix, if unspecified
    R   = get_matrix(ttT = ttT, which_matrix = which_matrix, which_tfbms = which_tfbms, explicit_zeros = explicit_zeros)
    ttT = append_matrix(ttT = ttT, ttT_sub = ttT_sub, R = R)

    return(ttT = ttT)
  }
