#' Get Power Estimates from EB functions
#'
#' @param pp_fun A power prior function from [binom.PP.EB()] with `mix=FALSE`)
#'
#' @return A `data.frame` containing the power prior weights
#' @export
#'
#' @examples
#' pp_eb <- binom.PP.EB(x = c(24,35), n = c(50,50), X = 0:50, N=50, mix = FALSE)
#' get_eb_powers(pp_eb)

get_eb_powers <- function(pp_fun) {
  stopifnot(is.function(pp_fun))
  e <- environment(pp_fun)
  if (!exists("ds", envir = e)) stop("Power prior weights are not defined.")
  deltas <- do.call(rbind.data.frame, e$ds)
  colnames(deltas) <- paste0("delta_", seq_len(NCOL(deltas)))
  cbind(X = e$X, deltas)
}
