library(sleuth)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]]
output <- args[[2]]

stab <- read.table(input, header=TRUE, stringsAsFactors=FALSE)
so <- sleuth_prep(stab)

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)
sleuth_significant <- dplyr::select(sleuth_significant, target_id, test_stat, pval, qval)

write.table(sleuth_significant, file=output, quote=FALSE, row.names=FALSE, append=TRUE, sep='\t')
