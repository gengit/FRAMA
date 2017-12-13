# checks GNU Parallel logs

args = commandArgs(trailingOnly = T)
log = read.table(args[1], header = T, sep = "\t")

if( nrow(log) > 0 & !is.null(log$Exitval) & sum(log$Exitval) == 0 ) {
  q("no", status =0)
} else {
  cat("Some jobs returned errors.")
  q("no", status = 1)
}

