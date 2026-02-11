setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library(Biostrings)
source("Functions.R")

# Perform the downstream analysis for differential expression analysis. (MaxLFQ)
analysis <- DownstreameR(report_file = "Data/20260202_085022_2041_Report_DownStreameR.tsv",
                         metadata_file = "Data/2041_Metadata.xlsx",
                         level = "protein", 
                         tool = "spectronaut",
                         quantification = "maxlfq",
                         group_replicates = TRUE,
                         protein_group_qvalue = 0.01,
                         remove_decoy = TRUE,
                         quantified_per_condition = 0,
                         perform_log2 = TRUE,
                         fasta_scale = TRUE,
                         perform_normalization = TRUE,
                         normalization_method = "median",
                         imputation = TRUE,
                         conditional_imputation = TRUE,
                         custom_condition_colors = c("#317EC2", "#C03830", "#008835"),
                         fdr = 0.05,
                         log2_fc = 0.58496250073,
                         report_only_regulated = FALSE)

# Get the protein expression df.
expressions <- analysis$Report$Data[ , c("protein_group", "protein_name", "fasta_file", analysis$Report$Quantitative_columns)]

# Map species to fasta files.
species_per_fasta <- setNames(nm = c("uniprotkb_proteome_UP000001360_Bifidobacterium_longum_2026_02_01",
                                     "uniprotkb_proteome_UP000008702_Bifidobacterium_adolescentis_2026_02_01",
                                     "uniprotkb_proteome_UP000008178_Roseburia_hominis_2026_02_01",
                                     "uniprotkb_proteome_UP000294398_Roseburia_intestinalis_2026_02_01",
                                     "uniprotkb_taxonomy_id_622312_Roseburia_inulinivorans_2026_02_01"),
                              object = c("Bifidobacterium longum subsp. infantis",
                                         "Bifidobacterium adolescentis",
                                         "Roseburia hominis",
                                         "Roseburia intestinalis",
                                         "Roseburia inulinivorans"))

# Keep only the species of major protein.
expressions$species <- unlist(x = lapply(expressions$fasta_file, function(x){
  
  files <- strsplit(x = x, split = ";")[[1]]
  
  species_per_file <- species_per_fasta[files]
  
  return(species_per_file[1])

}))

expressions <- expressions[ , c("protein_group", "protein_name", "species", analysis$Report$Quantitative_columns)]

##################################################################################
# Get how many proteins per species are in each sample.
##################################################################################
species_count_quantified_per_sample <- sapply(analysis$Report$Quantitative_columns, function(x){
  
  expr <- expressions[!is.na(x = expressions[ , x]),]
  
  species_count <- setNames(object = sapply(species_per_fasta, function(y){
    
    return(sum(expr$species == y, na.rm = TRUE))
  
  }), nm = species_per_fasta)
  
  return(species_count)
  
})

species_count_quantified_per_sample_long <- reshape2::melt(data = species_count_quantified_per_sample)
colnames(x = species_count_quantified_per_sample_long) <- c("Species", "Sample", "Count")

plot_species_count_bars <- ggplot2::ggplot(data = species_count_quantified_per_sample_long, mapping = ggplot2::aes(x = Sample, y = Count, fill = Species)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5) +
  ggplot2::geom_text(data = subset(x = species_count_quantified_per_sample_long, subset = Count != 0), mapping = ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.7), color = "white", size = 3.5, fontface = "bold") +
  ggplot2::scale_fill_manual(values = setNames(object = c("#FCB14F", "#F14040", "#7D2DB9", "#008000", "#1A6FDF"), nm = species_per_fasta)) +
  ggplot2::ylab("Protein Count") +
  ggplot2::xlab("Samples") +
  ggplot2::ggtitle("Protein count per species") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))

##################################################################################
# Find IBAQ proportions per species (stoichiometries).
##################################################################################
# Run partially the downstream analysis with iBAQ intensities.
ibaq_report <- import_report(file = "Data/20260202_085022_2041_Report_DownStreameR.tsv",
                             level = "protein", 
                             tool = "spectronaut",
                             quantification = "ibaq")

colnames(x = ibaq_report$Data)[8:16] <- analysis$Report$Quantitative_columns
ibaq_report$Quantitative_columns <- analysis$Report$Quantitative_columns

ibaq_report <- preprocessing(report = ibaq_report,
                             design = analysis$Metadata$Design,
                             group_replicates = TRUE,
                             protein_group_qvalue = 0.01,
                             remove_decoy = TRUE,
                             quantified_per_condition = 1,
                             perform_log2 = TRUE,
                             perform_normalization = TRUE,
                             normalization_method = "median",
                             imputation = FALSE,
                             conditional_imputation = FALSE)

ibaq_expressions <- ibaq_report$Data[ , c("protein_group", "protein_name", "fasta_file", ibaq_report$Quantitative_columns)]

# Keep only the species of major protein.
ibaq_expressions$species <- unlist(x = lapply(ibaq_expressions$fasta_file, function(x){
  
  files <- strsplit(x = x, split = ";")[[1]]
  
  species_per_file <- species_per_fasta[files]
  
  return(species_per_file[1])
  
}))

ibaq_expressions <- ibaq_expressions[ , c("protein_group", "protein_name", "species", ibaq_report$Quantitative_columns)]

ibaq_count_per_species_per_sample <- sapply(species_per_fasta, function(x){
  
  expr <- ibaq_expressions[ibaq_expressions$species == x, ibaq_report$Quantitative_columns]
  expr <- 2^expr
  
  return(colSums(x = expr, na.rm = TRUE))
  
})

colnames(x = ibaq_count_per_species_per_sample) <- species_per_fasta

# Un-log.
ibaq_totals <- colSums(x = 2^ibaq_expressions[ , ibaq_report$Quantitative_columns], na.rm = TRUE)

# Turn into percentages.
species_stoichiometry_per_sample <- t(100 *(ibaq_count_per_species_per_sample / ibaq_totals))

species_stoichiometry_per_sample_long <- reshape2::melt(data = species_stoichiometry_per_sample)
colnames(x = species_stoichiometry_per_sample_long) <- c("Species", "Sample", "Count")
species_stoichiometry_per_sample_long$Count <- round(x = species_stoichiometry_per_sample_long$Count, digits = 2) 

plot_stoichiometry_bars <- ggplot2::ggplot(data = species_stoichiometry_per_sample_long, mapping = ggplot2::aes(x = Sample, y = Count, fill = Species)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5) +
  ggplot2::geom_text(data = subset(x = species_stoichiometry_per_sample_long, subset = Count != 0), mapping = ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.7), color = "white", size = 3.5, fontface = "bold") +
  ggplot2::scale_fill_manual(values = setNames(object = c("#FCB14F", "#F14040", "#7D2DB9", "#008000", "#1A6FDF"), nm = species_per_fasta)) +
  ggplot2::ylab("Stoichiometry (%)") +
  ggplot2::xlab("Samples") +
  ggplot2::ggtitle("Species stoichiometry per sample") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))

# Export stoichiometries.
write.table(x = species_stoichiometry_per_sample_long, file = "Data/iBAQ_stoichiometries.csv", sep = ",", quote = FALSE)

##################################################################################
# Get how many proteins per species are expressed in each comparison.
##################################################################################
regulation_per_species_per_comparison <- do.call(what = "rbind", lapply(names(x = analysis$Stats$Stats), function(x){
  
  stats <- analysis$Stats$Stats[[x]]
  
  stats$Comparison <- rep(x = x, times = nrow(x = stats))
  
  return(stats[ , c("Comparison", "fasta_file", "Regulation")])
  
}))

#Keep only the species of major protein.
regulation_per_species_per_comparison$Species <- unlist(x = lapply(regulation_per_species_per_comparison$fasta_file, function(x){
  
  files <- strsplit(x = x, split = ";")[[1]]
  
  species_per_file <- species_per_fasta[files]
  
  return(species_per_file[1])
  
}))

# Split to down and up regulated bins.
regulation_per_species_per_comparison <- regulation_per_species_per_comparison[ , c("Comparison", "Regulation", "Species")]

downregulated_per_species_per_comparison <- regulation_per_species_per_comparison[regulation_per_species_per_comparison$Regulation == "Down-regulated", ]
rownames(x = downregulated_per_species_per_comparison) <- 1:nrow(x = downregulated_per_species_per_comparison)

upregulated_per_species_per_comparison <- regulation_per_species_per_comparison[regulation_per_species_per_comparison$Regulation == "Up-regulated", ]
rownames(x = upregulated_per_species_per_comparison) <- 1:nrow(x = upregulated_per_species_per_comparison)

# Barplot of downregulated protein count per species for each comparison.
downregulated_per_species_per_comparison_table <- t(table(downregulated_per_species_per_comparison$Comparison, downregulated_per_species_per_comparison$Species))
downregulated_per_species_per_comparison_table_long <- reshape2::melt(data = downregulated_per_species_per_comparison_table)
colnames(x = downregulated_per_species_per_comparison_table_long) <- c("Species", "Comparison", "Count")
downregulated_per_species_per_comparison_table_long$Comparison <- factor(downregulated_per_species_per_comparison_table_long$Comparison,
                                                                         levels = c("Mix_vs_HMOs", "Mix_vs_Fiber", "HMOs_vs_Fiber"))

plot_downregulated_bars <- ggplot2::ggplot(data = downregulated_per_species_per_comparison_table_long, mapping = ggplot2::aes(x = Comparison, y = Count, fill = Species)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5) +
  ggplot2::geom_text(data = subset(x = downregulated_per_species_per_comparison_table_long, subset = Count != 0), mapping = ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.7), color = "white", size = 3.5, fontface = "bold") +
  ggplot2::scale_fill_manual(values = setNames(object = c("#FCB14F", "#F14040", "#7D2DB9", "#008000", "#1A6FDF"), nm = species_per_fasta)) +
  ggplot2::ylab("Down-regulated Protein Count") +
  ggplot2::xlab("Comparison") +
  ggplot2::ggtitle("Number of down-regulated proteins per species for each comparison") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))

# Barplot of upregulated protein count per species for each comparison.
upregulated_per_species_per_comparison_table <- t(table(upregulated_per_species_per_comparison$Comparison, upregulated_per_species_per_comparison$Species))
upregulated_per_species_per_comparison_table_long <- reshape2::melt(data = upregulated_per_species_per_comparison_table)
colnames(x = upregulated_per_species_per_comparison_table_long) <- c("Species", "Comparison", "Count")
upregulated_per_species_per_comparison_table_long$Comparison <- factor(upregulated_per_species_per_comparison_table_long$Comparison,
                                                                         levels = c("Mix_vs_HMOs", "Mix_vs_Fiber", "HMOs_vs_Fiber"))

plot_upregulated_bars <- ggplot2::ggplot(data = upregulated_per_species_per_comparison_table_long, mapping = ggplot2::aes(x = Comparison, y = Count, fill = Species)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5) +
  ggplot2::geom_text(data = subset(x = upregulated_per_species_per_comparison_table_long, subset = Count != 0), mapping = ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.7), color = "white", size = 3.5, fontface = "bold") +
  ggplot2::scale_fill_manual(values = setNames(object = c("#FCB14F", "#F14040", "#7D2DB9", "#008000", "#1A6FDF"), nm = species_per_fasta)) +
  ggplot2::ylab("Up-regulated Protein Count") +
  ggplot2::xlab("Comparison") +
  ggplot2::ggtitle("Number of up-regulated proteins per species for each comparison") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))

#################################################################################
# Export the above figures.
#################################################################################
page1 <- (plot_species_count_bars / plot_stoichiometry_bars) + 
  patchwork::plot_annotation(title = "General species specific figures", 
                             caption = paste0("Page ", 1),
                             theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50)))

page2 <- (plot_downregulated_bars / plot_upregulated_bars) + 
  patchwork::plot_annotation(title = "General species specific figures", 
                             caption = paste0("Page ", 2),
                             theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50)))

pdf(file = "Data/General_figures.pdf", width = 8.27*1.4, height = 11.69*1.4)

invisible(x = print(page1))
invisible(x = print(page2))

dev.off()

##################################################################################
# Create fasta files for dbCAN3
##################################################################################
# Get and combine all reference proteomes.
fasta_paths <- paste0("Data/FASTA/", names(x = species_per_fasta), ".fasta")
fasta <- lapply(fasta_paths, function(x){return(Biostrings::readAAStringSet(filepath = x, format = "fasta"))})
fasta <- do.call(what = c, fasta)
fasta_names <- names(x = fasta)

# For comparison HMOs_vs_Fiber
HMOs_vs_Fiber <- analysis$Stats$Stats$HMOs_vs_Fiber
HMOs_vs_Fiber <- HMOs_vs_Fiber[HMOs_vs_Fiber$Regulation != "Non-regulated", ]
rownames(x = HMOs_vs_Fiber) <- 1:nrow(x = HMOs_vs_Fiber)
HMOs_vs_Fiber$protein_group <- sapply(HMOs_vs_Fiber$protein_group, function(x){return(strsplit(x = x, split = ";")[[1]][1])})

# Get and save the sequences of the regulated proteins.
regulated_idx <- sapply(unlist(x = HMOs_vs_Fiber$protein_group), function(y){return(which(x = grepl(pattern = y, x = fasta_names)))})
HMOs_vs_Fiber_fasta <- fasta[regulated_idx]
names(x = HMOs_vs_Fiber_fasta) <- gsub(pattern = "->", replacement = "-", x = names(x = HMOs_vs_Fiber_fasta))
Biostrings::writeXStringSet(HMOs_vs_Fiber_fasta, paste0("Data\\dbCAN3\\HMOs_vs_Fiber.fasta"), format = "fasta", width = 200001)

# For comparison Mix_vs_Fiber
Mix_vs_Fiber <- analysis$Stats$Stats$Mix_vs_Fiber
Mix_vs_Fiber <- Mix_vs_Fiber[Mix_vs_Fiber$Regulation != "Non-regulated", ]
rownames(x = Mix_vs_Fiber) <- 1:nrow(x = Mix_vs_Fiber)
Mix_vs_Fiber$protein_group <- sapply(Mix_vs_Fiber$protein_group, function(x){return(strsplit(x = x, split = ";")[[1]][1])})

# Get and save the sequences of the regulated proteins.
regulated_idx <- sapply(unlist(x = Mix_vs_Fiber$protein_group), function(y){return(which(x = grepl(pattern = y, x = fasta_names)))})
Mix_vs_Fiber_fasta <- fasta[regulated_idx]
names(x = Mix_vs_Fiber_fasta) <- gsub(pattern = "->", replacement = "-", x = names(x = Mix_vs_Fiber_fasta))
Biostrings::writeXStringSet(Mix_vs_Fiber_fasta, paste0("Data\\dbCAN3\\Mix_vs_Fiber.fasta"), format = "fasta", width = 200001)

# For comparison Mix_vs_HMOs
Mix_vs_HMOs <- analysis$Stats$Stats$Mix_vs_HMOs
Mix_vs_HMOs <- Mix_vs_HMOs[Mix_vs_HMOs$Regulation != "Non-regulated", ]
rownames(x = Mix_vs_HMOs) <- 1:nrow(x = Mix_vs_HMOs)
Mix_vs_HMOs$protein_group <- sapply(Mix_vs_HMOs$protein_group, function(x){return(strsplit(x = x, split = ";")[[1]][1])})

# Get and save the sequences of the regulated proteins.
regulated_idx <- sapply(unlist(x = Mix_vs_HMOs$protein_group), function(y){return(which(x = grepl(pattern = y, x = fasta_names)))})
Mix_vs_HMOs_fasta <- fasta[regulated_idx]
names(x = Mix_vs_HMOs_fasta) <- gsub(pattern = "->", replacement = "-", x = names(x = Mix_vs_HMOs_fasta))
Biostrings::writeXStringSet(Mix_vs_HMOs_fasta, paste0("Data\\dbCAN3\\Mix_vs_HMOs.fasta"), format = "fasta", width = 200001)

