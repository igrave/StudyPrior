get_powers <- function(pp_fun) {
   e <- environment(pp_fun)
   deltas <- do.call(rbind.data.frame, e$ds)
   colnames(deltas) <- paste0("delta_", seq_len(NCOL(deltas)))
   cbind(X = e$X, deltas)
}
