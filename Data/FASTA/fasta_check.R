setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(Biostrings)

fasta_622312 <- Biostrings::readAAStringSet(filepath = "uniprotkb_taxonomy_id_622312_Roseburia_inulinivorans_2026_02_01.fasta", format = "fasta")
proteins_622312 <- names(x = fasta_622312)

fasta_622312_2 <- Biostrings::readAAStringSet(filepath = "uniprotkb_proteome_UP000003561_2026_02_02.fasta", format = "fasta")
proteins_622312_2 <- names(x = fasta_622312_2)


setdiff(x = proteins_622312, y = proteins_622312_2)

