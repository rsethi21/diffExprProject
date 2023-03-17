# loading libraries assuming already installed locally
library(sleuth)
library(dplyr)

# args input form the command line to allow for flexibility of code
args <- commandArgs(trailingOnly = TRUE)
input <- args[[1]] # input of TPMs
output <- args[[2]] # output path
intermediate <- args[[3]] # intermediate file that saves this output for a different program

stab <- read.table(input, header=TRUE, stringsAsFactors=FALSE) # sleuth table to read in TPMs
so <- sleuth_prep(stab) # create sleuth object

# fit sleuth object to run statistical analysis and find differentially expressed genes between samples using the TPMs calculated from kallisto
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) # finding significantly different genes expressed betweeen samples of FDR below 0.05
sleuth_significant <- dplyr::select(sleuth_significant, target_id, test_stat, pval, qval) # saving the id, test statistic, p and q values of the differentially expressed genes found most significant

# write the output to the specified folder
write.table(sleuth_significant, file=output, quote=FALSE, row.names=FALSE, append=TRUE, sep='\t')
write.table(sleuth_significant, file=intermediate, quote=FALSE, row.names=FALSE, append=TRUE, sep='\t')
