setwd()
set.seed(42)

library(dplyr)

# import raw
raw_counts_mat = read.delim("../99.processed/counts_mat.txt", sep = '\t', check.names = F)

# cleaning
counts_mat = as.matrix(counts_mat)

# extract metadata
meta = data.frame(orig.ident = colnames(counts_mat), stringsAsFactors = F)
meta$CellLine = gsub("MgEtOH", "", meta$orig.ident)
meta$concentration = gsub(".*\\_", "", meta$orig.ident)
meta$PRS = ifelse(meta$CellLine %in% c("100", "203", "204", "528", "864"), "High PRS", "Low PRS")

# write out
write.table(counts_mat, "XXX/counts_mat.txt", quote = F)
write.table(meta, "XXX/meta.txt", quote = F, row.names = F)