# Dependencies
packages <- c("checkmate", "readxl", "dplyr", "tidyr", "openxlsx", "reshape2", "ggplot2", "gridExtra", "patchwork")
bioconductor <- c("limma")

# Install missing dependencies.
to_install <- setdiff(x = packages, y = rownames(x = installed.packages()))
to_install_b <- setdiff(x = bioconductor, y = rownames(x = installed.packages()))

if(length(x = to_install) > 0){
  
  install.packages(to_install)
  
}

if(length(x = to_install_b) > 0){
  
  if (!require("BiocManager", quietly = TRUE)){
    
    install.packages("BiocManager")
    
  }
  
  BiocManager::install(to_install_b)
  
}

# Load libraries
library(checkmate)
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)
library(limma)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(patchwork)

DownstreameR <- function(report_file = NULL, 
                         metadata_file = NULL, 
                         level = "protein", 
                         tool = "spectronaut", 
                         quantification = "maxlfq",
                         group_replicates = TRUE,
                         protein_group_qvalue = 0.01,
                         remove_decoy = TRUE,
                         quantified_per_condition = 0,
                         perform_log2 = TRUE,
                         fasta_scale = FALSE,
                         perform_normalization = TRUE,
                         normalization_method = "median",
                         imputation = FALSE,
                         conditional_imputation = FALSE,
                         custom_condition_colors = FALSE,
                         fdr = 0.05,
                         log2_fc = 0.58496250073,
                         report_only_regulated = TRUE){
  
  # Validate fdr (must be a single numeric between 0 and 1).
  validate_fdr <- checkmate::check_number(x = fdr, lower = 0, upper = 1, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_fdr)){
    
    fdr <- 0.05
    warning("Warning at function DownstreameR():\n\tInvalid fdr value. Default of 0.05 is set.")
    
  }
  
  # Validate log2_fc (must be a single positive numeric).
  validate_log2_fc <- checkmate::check_number(x = log2_fc, lower = 0, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_log2_fc)){
    
    log2_fc <- 0.58496250073
    warning("Warning at function DownstreameR():\n\tInvalid log2_fc value. Default of 0.58496250073 is set.")
    
  }
  
  # Validate report_only_regulated argument.
  validate_report_only_regulated <- checkmate::check_flag(x = report_only_regulated, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_report_only_regulated)){
    
    group_replicates <- TRUE
    warning("Warning at function DownstreameR():\n\tInvalid report_only_regulated value. TRUE is set as default.")
    
  }
  
  # Load the report file.
  cat("1. Importing report file\n")
  report <- import_report(file = report_file, level = level, tool = tool, quantification = quantification)
  if(is.null(x = report)){
    
    cat(paste0("   File ", report_file, " is either empty or level and tool options are not supported yet.\n\n"))
    return(NULL)
    
  }
  cat(paste0("   File ", report_file, " is succesfully imported.\n\n"))
  
  # Load the metadata file.
  cat("2. Importing metadata file\n")
  metadata <- import_metadata(file = metadata_file)
  cat(paste0("   File ", metadata_file, " is succesfully imported.\n\n"))
  
  # Add sample names to the corresponding quantitative columns in data.
  cat("3. Aligning metadata to experimental data\n")
  sample_column_indices <- sapply(metadata$Samples$`Raw File`, function(x){
    
    index = which(x = grepl(pattern = x, x = colnames(x = report$Data)))
    if(length(x = index) == 0){
      
      index <- NA
      
    }
    
    return(index)
    
  })
  
  # If no quantitative columns are found for all raw files raise an error.
  # Otherwise rename the quantitative data columns according to sample names and keep only the quantitative columns corresponding to a sample found in metadata.
  if(any(is.na(x = sample_column_indices))){
    
    stop(paste0("Error at function DownstreameR():\n\tNo quantitative columns found in file ", report_file, " corresponding to raw file(s) ",
                paste(names(x = sample_column_indices)[is.na(x = sample_column_indices)], collapse = ", "), "."))
    
  } else {
    
    colnames(x = report$Data)[sample_column_indices] <- metadata$Samples$`Sample Name`
    
    report$Data <- report$Data[ , c(report$Standard_columns, metadata$Samples$`Sample Name`)]
    
    quantitative_column_difference <- length(x = report$Quantitative_columns) - length(x = metadata$Samples$`Sample Name`)
    
    if(quantitative_column_difference > 0){
      
      cat(paste0("   Omitted ", quantitative_column_difference, " quantitative columns from the report, since they couldn't be mapped to the metadata raw files.\n"))
      
    }
    
    report$Quantitative_columns <- metadata$Samples$`Sample Name`
    
  }
  
  cat("   Quantitative data found for all raw files.\n\n")
  
  # Preproccessing of report.
  cat("4. Filtering and preprocessing report\n")
  report <- preprocessing(report = report,
                          design = metadata$Design,
                          group_replicates = group_replicates,
                          protein_group_qvalue = protein_group_qvalue,
                          remove_decoy = remove_decoy,
                          quantified_per_condition = quantified_per_condition,
                          perform_log2 = perform_log2,
                          fasta_scale = fasta_scale,
                          perform_normalization = perform_normalization,
                          normalization_method = normalization_method,
                          imputation = imputation,
                          conditional_imputation = conditional_imputation)
  
  cat("   Preprocessing is successfully done.\n")
  
  # Save the processed report in xlsx format.
  processed_report_file <- paste0(dirname(path = report_file), "/DownstreameR_report.xlsx")
  
  # Create workbook and add a worksheet
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Data")
  headerStyle <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = c("top", "bottom", "left", "right"))
  openxlsx::writeData(wb = wb, sheet = "Data", x = report$Data, headerStyle = headerStyle)
  openxlsx::setColWidths(wb = wb, sheet = "Data", cols = 1:ncol(x = report$Data), widths = "auto")
  
  openxlsx::addWorksheet(wb = wb, sheetName = "Samples_metadata")
  openxlsx::writeData(wb = wb, sheet = "Samples_metadata", x = metadata$Samples, headerStyle = headerStyle)
  openxlsx::setColWidths(wb = wb, sheet = "Samples_metadata", cols = 1:ncol(x = metadata$Samples), widths = "auto")
  
  openxlsx::saveWorkbook(wb, processed_report_file, overwrite = TRUE)
  
  cat(paste0("   Metadata and preprocessed data are exported at: ", processed_report_file, ".\n"))
  
  # Create QC plots
  qc_plots <- protein_qc_plots(report = report,
                               metadata = metadata,
                               perform_log2 = perform_log2,
                               normalization_method = normalization_method,
                               custom_condition_colors = custom_condition_colors)
  
  report[["QC_plots"]] <- qc_plots
  
  cat("   QC plots are generated.\n\n")
  
  # Statistics.
  cat("5. Statistical analysis\n")
  if(metadata$Do_statistics){
    
    statistics <- perform_statistics(report = report, metadata = metadata)
    
    if(is.null(x = statistics)){
      
      cat("   No statistics are calculated. Check the warning message above to find the reason.\n\n")
      
    } else {
      
      cat("   Statistical analysis is completed.\n")
      
      # Filter for fdr and log2_fc and prepare for visualization.
      statistics$Stats <- lapply(statistics$Stats, function(x){
        
        filter_down <- (x$FDR <= fdr) & (x$logFC <= -log2_fc)
        filter_down[is.na(x = filter_down)] <- FALSE
        filter_up <- (x$FDR <= fdr) & (x$logFC >= log2_fc)
        filter_up[is.na(x = filter_up)] <- FALSE
        
        x$Regulation <- rep(x = NA, times = nrow(x = x))
        x$Regulation[filter_down] <- "Down-regulated"
        x$Regulation[filter_up] <- "Up-regulated"
        x$Regulation[is.na(x = x$Regulation)] <- "Non-regulated"
        
        x$Keep <- filter_down | filter_up
        
        return(x)
      
      })
      
      # Create figures per comparison.
      statistics[["Statistics_plots"]] <- plot_statistics(statistics = statistics$Stats, fdr = fdr, log2_fc = log2_fc)
      
      # Add the protein information and prepare for report.
      statistics$Stats <- lapply(statistics$Stats, function(x){

        x <- cbind(report$Data[ , report$Standard_columns], x)
        
        if(report_only_regulated){
          
          x <- x[x$Keep , ]
        
        }
        
        x <- x[ , -which(x = colnames(x = x) == "Keep")]

        if(nrow(x = x) > 0){

          rownames(x = x) <- 1:nrow(x = x)

        }

        return(x)

      })
      
      # Add the contrasts and statistics in the report. A sheet per comparison.
      wb <- openxlsx::loadWorkbook(file = processed_report_file)
      
      openxlsx::addWorksheet(wb = wb, sheetName = "Contrasts")
      openxlsx::writeData(wb = wb,
                          sheet = "Contrasts",
                          x = data.frame(Conditions = rownames(x = statistics$Contrasts), statistics$Contrasts, stringsAsFactors = FALSE),
                          headerStyle = headerStyle)
      openxlsx::setColWidths(wb = wb, sheet = "Contrasts", cols = 1:(ncol(x = statistics$Contrasts) + 1), widths = "auto")
      
      for(comparison in names(x = statistics$Stats)){
        
        openxlsx::addWorksheet(wb = wb, sheetName = comparison)
        openxlsx::writeData(wb = wb, sheet = comparison, x = statistics$Stats[[comparison]], headerStyle = headerStyle)
        openxlsx::setColWidths(wb = wb, sheet = comparison, cols = 1:ncol(x = statistics$Stats[[comparison]]), widths = "auto")
        
      }
      
      openxlsx::saveWorkbook(wb, processed_report_file, overwrite = TRUE)
      
      cat(paste0("   A sheet per comparison is added in file ", processed_report_file, ".\n\n"))
      
    }
    
  } else {
    
    cat("   No statistics are calculated, None column is used in the metadata contrast table.\n\n")
    statistics <- NULL
    
  }
  
  cat("6. Processing graphics\n")
  # Create the QC pages of the report.
  report_pages <- list(QC_1 = (report$QC_plots$NA_plot / report$QC_plots$Distribution_plot) + 
                         patchwork::plot_annotation(title = "Quality control",
                                                    subtitle = "Missing value and distributions",
                                                    caption = "Page 1",
                                                    theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50))),
                       QC_2 = (report$QC_plots$Violin_plot / report$QC_plots$CV_plot) +
                         patchwork::plot_annotation(title = "Quality control",
                                                    subtitle = "Reproducibility assessment",
                                                    caption = "Page 2",
                                                    theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50))),
                       QC_3 = (report$QC_plots$Heatmap / report$QC_plots$PCA_plot) +
                         patchwork::plot_annotation(title = "Quality control",
                                                    subtitle = "Sample-wise similarity assessment",
                                                    caption = "Page 3",
                                                    theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50))))
  # Add the pages of the statistics figures.
  if(!is.null(x = statistics)){
    
    report_pages <- c(report_pages, statistics$Statistics_plots)
    
  }
  
  cat("   Report pages are generated.\n")
  
  # Create file name.
  figure_report_file <- paste0(dirname(path = report_file), "/DownstreameR_figures.pdf")
  pdf(file = figure_report_file, width = 8.27*1.4, height = 11.69*1.4)
  
  for(page in report_pages) {
    
    invisible(x = print(page))
    
  }
  
  dev.off()
  
  cat(paste0("   QC and statistical analysis figures are added in: ", figure_report_file, ".\n\n"))
  
  # Get function arguments and their values, and add at the report.
  cat("7. Saving parameters\n")
  arg_names <- names(x = formals())
  arg_values <- mget(x = arg_names, envir = environment())
  
  arg_table <- data.frame(Parameter = names(x = arg_values), Value = sapply(arg_values, function(x){ return(paste0(x, collapse = ", "))}), stringsAsFactors = FALSE)
  
  wb <- openxlsx::loadWorkbook(file = processed_report_file)
  openxlsx::addWorksheet(wb = wb, sheetName = "Parameters")
  openxlsx::writeData(wb = wb, sheet = "Parameters",x = arg_table, headerStyle = headerStyle)
  openxlsx::setColWidths(wb = wb, sheet = "Parameters", cols = 1:2, widths = "auto")
  openxlsx::saveWorkbook(wb, processed_report_file, overwrite = TRUE)
  
  cat(paste0("   Parameters are exported at: ", processed_report_file, ".\n\n"))
  
  cat("DONE!")

  return(list(Report = report, Metadata = metadata, Stats = statistics))
  
}

import_report <- function(file = NULL, level = "protein", tool = "spectronaut", quantification = "maxlfq"){
  
  argument_validation <- checkmate::makeAssertCollection()
  
  # Validate file argument.
  checkmate::assert(
    checkmate::check_string(x = file, na.ok = FALSE, null.ok = FALSE),
    checkmate::check_file_exists(x = file),
    combine = "and",
    add = argument_validation
  )
  
  # Validate level argument.
  checkmate::assert(
    checkmate::check_string(x = level, na.ok = FALSE, null.ok = FALSE),
    checkmate::check_choice(x = level, choices = c("protein", "peptide")),
    combine = "and",
    add = argument_validation
  )
  
  # Validate tool argument.
  checkmate::assert(
    checkmate::check_string(x = tool, na.ok = FALSE, null.ok = FALSE),
    checkmate::check_choice(x = tool, choices = c("spectronaut", "pd")),
    combine = "and",
    add = argument_validation
  )
  
  # Validate quantification argument.
  checkmate::assert(
    checkmate::check_string(x = quantification, na.ok = FALSE, null.ok = FALSE),
    checkmate::check_choice(x = quantification, choices = c("maxlfq", "ibaq")),
    combine = "and",
    add = argument_validation
  )
  
  if(!argument_validation$isEmpty()){
    
    stop(paste(c("Error at function import_report():", paste0("\t", argument_validation$getMessages())), collapse = "\n"))
    
  }
  
  # Define mandatory columns.
  spectronaut_fixed_protein_columns <- c("PG.Qvalue",
                                         "PG.ProteinGroups",
                                         "PG.ProteinAccessions",
                                         "PG.ProteinDescriptions",
                                         "PG.ProteinNames",
                                         "PG.FastaFiles",
                                         "PG.NrOfPrecursorsIdentified..Experiment.wide.")
  
  spectronaut_pattern_columns <- ".PG.Quantity$|.PG.IBAQ$"
  
  spectronaut_fixed_protein_columns_rnm <- c("qvalue",
                                             "protein_group",
                                             "protein_accession",
                                             "protein_description", 
                                             "protein_name",
                                             "fasta_file",
                                             "unique_precursors")
  
  if(tool == "spectronaut"){
    
    if(level == "protein"){
      
      # Import data.
      data <- read.delim2(file = file,
                          header = TRUE,
                          na.strings = c("", "NA", "NaN"),
                          sep = "\t",
                          dec = ",",
                          stringsAsFactors = FALSE)
      
      # Check if there are data.
      if(nrow(x = data) > 0){
        
        # Check for the mandatory columns.
        columns_not_found <- setdiff(x = spectronaut_fixed_protein_columns, y = colnames(x = data))
        
        # If all mandatory columns are there check the suffixes of the quantitative columns.
        if(length(x = columns_not_found) == 0){
          
          columns_for_pattern_check <- setdiff(x = colnames(x = data), y = spectronaut_fixed_protein_columns)
          pattern_check <- grepl(pattern = spectronaut_pattern_columns, x = columns_for_pattern_check)
          
          # If all quantitative columns are according to the suffix patterns, keep only those corresponding to the quantification argument.
          if(all(pattern_check, na.rm = TRUE)){
            
            if(quantification == "maxlfq"){
              
              quantitative_column_mask <- grepl(pattern = strsplit(x = spectronaut_pattern_columns, split = "\\|")[[1]][1], x = columns_for_pattern_check)
              
            } else {
              
              quantitative_column_mask <- grepl(pattern = strsplit(x = spectronaut_pattern_columns, split = "\\|")[[1]][2], x = columns_for_pattern_check)
              
            }
            
            quantitative_columns <- columns_for_pattern_check[quantitative_column_mask]
            
            # Keep only mandatory and quantitative columns.
            data <- data[ , c(spectronaut_fixed_protein_columns, quantitative_columns)]
            
            # Make sure columns are in the right data type.
            data[ , spectronaut_fixed_protein_columns[1]] <- as.numeric(x = data[ , spectronaut_fixed_protein_columns[1]])
            data[ , spectronaut_fixed_protein_columns[-c(1, length(x = spectronaut_fixed_protein_columns))]] <- apply(data[ , spectronaut_fixed_protein_columns[-c(1, length(x = spectronaut_fixed_protein_columns))]], 2, as.character)
            data[ , spectronaut_fixed_protein_columns[length(x = spectronaut_fixed_protein_columns)]] <- as.integer(x = data[ , spectronaut_fixed_protein_columns[length(x = spectronaut_fixed_protein_columns)]])
            
            # Check if quantitative columns are numeric.
            # In case of ibaq, first keep only the major group ibaq scores and then cast to numeric.
            # In case of maxlfq, just cast to numeric.
            character_column_mask <- apply(data[ , quantitative_columns], 2, class) == "character"

            if(any(character_column_mask, na.rm = TRUE)){
              
              if(quantification == "ibaq"){
                
                data[ , quantitative_columns[character_column_mask]] <- apply(data[ , quantitative_columns[character_column_mask]], 2, function(x){ return(sapply(strsplit(x = x, split = ";"), function(x){ return(x[1])}))})
              
              }
                
              data[ , quantitative_columns[character_column_mask]] <- apply(data[ , quantitative_columns[character_column_mask]], 2, function(x){ return(as.numeric(x = gsub(pattern = ",", replacement = ".", x = x)))})
              
            }              
            
            # Rename columns.
            colnames(x = data) <- c(spectronaut_fixed_protein_columns_rnm,
                                    gsub(pattern = spectronaut_pattern_columns, replacement = "", x = quantitative_columns))
            
            report <- list(Data = data, Standard_columns = spectronaut_fixed_protein_columns_rnm, Quantitative_columns = quantitative_columns)
            
          } else {
            
            stop(paste0("Error at function import_report():\n\tThe quantitative column(s) ",
                        paste(columns_for_pattern_check[!pattern_check], collapse = ", "),
                        " do not have the .PG.Quantity or .PG.IBAQ suffix."))
            
          }
          
        } else {
          
          stop(paste0("Error at function import_report():\n\tMandatory column(s) ",
                      paste(columns_not_found, collapse = ", "), " not found in file ", file, "."))
          
        }
        
      } else {
        
        warning(paste0("Warning at function import_report():\n\tFile ", file, " is empty, a NULL object is returned."))
        report <- NULL
        
      }
      
    } else {
      
      # Placeholder: Peptide option to be added.
      warning("Warning at function import_report():\n\tPeptide reports from Spectronaut are not supported yet. A NULL object is returned.")
      report <- NULL
      
    }
    
  } else {
    
    # Placeholder: PD option to be added.
    warning("Warning at function import_report():\n\tPD reports are not supported yet. A NULL object is returned.")
    report <- NULL
    
  }
  
  return(report)
  
}

import_metadata <- function(file = NULL){
  
  # Validate file argument.
  file_validation <- checkmate::makeAssertCollection()
  
  checkmate::assert(
    checkmate::check_string(x = file, na.ok = FALSE, null.ok = FALSE),
    checkmate::check_file_exists(x = file, extension = ".xlsx"),
    combine = "and",
    add = file_validation
  )
  
  if(!file_validation$isEmpty()){
    
    stop(paste(c("Error at function import_metadata():", paste0("\t", file_validation$getMessages())), collapse = "\n"))
    
  }
  
  # Define constrains.
  na_strings <- c("", " ", "NaN", "NA")
  metadata_sheet_names <- c("Samples", "Contrasts")
  mandatory_samples_sheet_columns <- c("Sample ID", "Sample Name", "Raw File", "Condition", "Replicate")
  
  # Check if the imported file can be loaded.
  error <- suppressWarnings(expr = try(expr = as.data.frame(x = readxl::read_excel(path = file, sheet = 1, .name_repair = "minimal"), stringsAsFactors = FALSE), silent = TRUE))
  
  if(class(x = error) == "try-error"){
  
    stop(paste0("Error at function import_metadata():\n\tThe file ", file, " could not be loaded."))      
    
  } else {
    
    # Read the sheet names.
    sheet_names <- readxl::excel_sheets(path = file)
    
    # Check if there are both the "Samples" and "Contrast" sheets.
    missing_sheets <- setdiff(x = metadata_sheet_names, y = sheet_names)
    
    if(length(x = missing_sheets) == 0){
      
      # Read the columns of the Samples sheet.
      samples_sheet_columns <- names(x = readxl::read_excel(path = file, sheet = metadata_sheet_names[1], .name_repair = "minimal", n_max = 0))
      
      # Validate Samples sheet columns.
      samples_sheet_columns_not_found <- setdiff(x = mandatory_samples_sheet_columns, y = samples_sheet_columns)
      
      if(length(x = samples_sheet_columns_not_found) == 0){
        
        # Import Samples sheet.
        samples_sheet <- as.data.frame(x = readxl::read_excel(path = file,
                                                              sheet = metadata_sheet_names[1],
                                                              range = readxl::cell_cols(x = 1:length(x = samples_sheet_columns)),
                                                              col_types = c("numeric", "text", "text", "text", "numeric"),
                                                              na = na_strings, 
                                                              .name_repair = "minimal"), stringsAsFactors = FALSE)
        
        # There should not be any missing values.
        samples_number_of_nas <- apply(samples_sheet, 2, function(x){ return(sum(is.na(x = x), na.rm = TRUE))})
        
        if(any(samples_number_of_nas > 0, na.rm = TRUE)){
          
          stop(paste0("Error at function import_metadata():\n\tThe column(s) ", paste(names(x = samples_number_of_nas[samples_number_of_nas > 0]), collapse = ", "),
                      " of Samples sheet contain missing values."))
          
        }
        
        # Create an experimental setup data frame of sample names per condition.
        experimental_design <- samples_sheet %>% 
          dplyr::select(Condition, Replicate, "Sample Name") %>% 
          tidyr::pivot_wider(names_from = Condition, values_from = "Sample Name", values_fill = NA) %>%
          dplyr::arrange(Replicate) %>%
          dplyr::select(-Replicate)
        
      } else {
        
        stop(paste0("Error at function import_metadata():\n\tThe Samples sheet does not contain the mandatory column(s) ",
                    paste(samples_sheet_columns_not_found, collapse = ", "), ".")) 
        
      }
      
      # Read the columns for the Contrasts sheet.
      contrast_sheet_columns <- names(x = readxl::read_excel(path = file, sheet = metadata_sheet_names[2], .name_repair = "minimal", n_max = 0))
      
      if(contrast_sheet_columns[1] == "Conditions"){
        
        if(contrast_sheet_columns[2] == "All"){
          
          if(length(x = contrast_sheet_columns) > 2){
            
            stop("Error at function import_metadata():\n\tThe Contrasts sheet when All column is used may have only 2 columns (Conditions and All).")
            
          }
          
          do_statistics <- TRUE
          
        } else if (contrast_sheet_columns[2] == "None"){
          
          if(length(x = contrast_sheet_columns) > 2){
            
            stop("Error at function import_metadata():\n\tThe Contrasts sheet when None column is used may have only 2 columns (Conditions and None).")
            
          }
          
          do_statistics <- FALSE
          
        } else {
          
          do_statistics <- TRUE
          
        }
        
        # Import Contrasts sheet.
        contrasts_sheet <- as.data.frame(x = readxl::read_excel(path = file,
                                                                sheet = metadata_sheet_names[2],
                                                                range = readxl::cell_cols(x = 1:length(x = contrast_sheet_columns)),
                                                                col_types = c("text", rep(x = "numeric", times = length(x = contrast_sheet_columns) - 1)),
                                                                na = na_strings, 
                                                                .name_repair = "minimal"), stringsAsFactors = FALSE)
        
        #Validate values of condition contrasts. There should not be NA and only the values 0, 1, and 2 are allowed.
        contrasts_invalid_values <- apply(contrasts_sheet, 2, function(x){
          
          return(c(na_values = sum(is.na(x = x), na.rm = TRUE),
                   invalid_values = sum(!(x %in% c(0, 1, 2)), na.rm = TRUE)))
          
        })
        
        # Validate NA values.
        if(any(contrasts_invalid_values["na_values", ] > 0, na.rm = TRUE)){
          
          stop(paste0("Error at function import_metadata():\n\tThe column(s) ", paste(colnames(x = contrasts_invalid_values)[contrasts_invalid_values["na_values", ] > 0], collapse = ", "),
                      " of Contrasts sheet contain missing values."))
          
        }
        
        # Validate contrasts values.
        if(any(contrasts_invalid_values["invalid_values", 2:ncol(x = contrasts_invalid_values)] > 0, na.rm = TRUE)){
          
          stop(paste0("Error at function import_metadata():\n\tThe column(s) ", 
                      paste(colnames(x = contrasts_invalid_values[ , 2:ncol(x = contrasts_invalid_values)])[contrasts_invalid_values["invalid_values", 2:ncol(x = contrasts_invalid_values)] > 0], collapse = ", "),
                      " of Contrasts sheet contain invalid values. Values must be 0, 1 or 2."))
          
        }
        
      } else {
        
        stop("Error at function import_metadata():\n\tThe Contrasts sheet does not contain the mandatory column Conditions.")
        
      }
      
      return(list(Samples = samples_sheet, Design = experimental_design, Contrasts = contrasts_sheet, Do_statistics = do_statistics))
      
    } else {
      
      stop(paste0("Error at function import_metadata():\n\tThe file ", file, " does not contain the mandatory sheet(s) ",
                  paste(missing_sheets, collapse = ", "), "."))   
      
    }
    
  }
  
}

preprocessing <- function(report = NULL,
                          design = NULL,
                          group_replicates = TRUE,
                          protein_group_qvalue = 0.01,
                          remove_decoy = TRUE,
                          quantified_per_condition = 0,
                          perform_log2 = TRUE,
                          fasta_scale = FALSE,
                          perform_normalization = TRUE,
                          normalization_method = "median",
                          imputation = FALSE,
                          conditional_imputation = FALSE){
  
  # Validate group_replicates argument.
  validate_group_replicates <- checkmate::check_flag(x = group_replicates, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_group_replicates)){
    
    group_replicates <- TRUE
    warning("Warning at function preprocessing():\n\tInvalid group_replicates value. TRUE is set as default.")
    
  }
  
  # Validate protein_group_qvalue (must be a single numeric between 0 and 1).
  validate_protein_group_qvalue <- checkmate::check_number(x = protein_group_qvalue, lower = 0, upper = 1, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_protein_group_qvalue)){
    
    protein_group_qvalue <- 0.01
    warning("Warning at function preprocessing():\n\tInvalid protein_group_qvalue value. Default of 0.01 is set.")
    
  }
  
  # Validate remove_decoy argument.
  validate_remove_decoy <- checkmate::check_flag(x = remove_decoy, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_remove_decoy)){
    
    remove_decoy <- TRUE
    warning("Warning at function preprocessing():\n\tInvalid remove_decoy value. TRUE is set as default.")
    
  }
  
  # Validate quantified_per_condition argument (must be a positive integer).
  validate_qauntified_per_condition <- check_int(x = quantified_per_condition, lower = 0, na.ok = FALSE)
  
  if(!isTRUE(x = validate_qauntified_per_condition)){
    
    quantified_per_condition <- 0
    warning("Warning at function preprocessing():\n\tInvalid quantified_per_condition value. Default of 0 is set.")
    
  }
  
  # Validate fasta_scale argument.
  validate_fasta_scale <- checkmate::check_flag(x = fasta_scale, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_fasta_scale)){
    
    perform_log2 <- FALSE
    warning("Warning at function preprocessing():\n\tInvalid fasta_scale value. FALSE is set as default.")
    
  }
  
  # Validate perform_log2 argument.
  validate_perform_log2 <- checkmate::check_flag(x = perform_log2, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_perform_log2)){
    
    perform_log2 <- TRUE
    warning("Warning at function preprocessing():\n\tInvalid perform_log2 value. TRUE is set as default.")
    
  }
  
  # Validate perform_normalization argument.
  validate_perform_normalization <- checkmate::check_flag(x = perform_normalization, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_perform_normalization)){
    
    perform_normalization <- TRUE
    warning("Warning at function preprocessing():\n\tInvalid perform_normalization value. TRUE is set as default.")
    
  }
  
  #Validate normalization_method argument.
  validate_normalization_method <- checkmate::makeAssertCollection()
  
  checkmate::assert(
    checkmate::check_string(x = normalization_method, na.ok = FALSE, null.ok = FALSE),
    checkmate::check_choice(x = normalization_method, choices = c("median", "mean")),
    combine = "and",
    add = validate_normalization_method
  )
  
  if(!validate_normalization_method$isEmpty()){
    
    normalization_method <- "median"
    warning("Warning at function preprocessing():\n\tInvalid normalization_method value. Default median normalization is set.")
    
  }
  
  # Validate imputation argument.
  validate_imputation <- checkmate::check_flag(x = imputation, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_imputation)){
    
    imputation <- FALSE
    warning("Warning at function preprocessing():\n\tInvalid imputation value. FALSE is set as default.")
    
  }
  
  # Validate conditional imputation argument.
  validate_conditional_imputation <- checkmate::check_flag(x = conditional_imputation, na.ok = FALSE, null.ok = FALSE)
  
  if(!isTRUE(x = validate_conditional_imputation)){
    
    conditional_imputation <- FALSE
    warning("Warning at function preprocessing():\n\tInvalid conditional_imputation value. FALSE is set as default.")
    
  }
  
  # Create a vector of quantitative column names ordered by condition.
  ordered_quantitative_columns <- na.omit(object = unlist(x = design, use.names = FALSE))
  
  # Q-value filtering.
  report$Data <- report$Data[report$Data$qvalue <= protein_group_qvalue, ]
  
  # Remove decoys.
  report$Data <- report$Data[!grepl(pattern = "_Decoy$", x = report$Data$protein_group), ]
  
  # Arrange columns based on conditions.
  if(group_replicates){
    
    report$Data <- report$Data[ , c(report$Standard_columns, ordered_quantitative_columns)]
    report$Quantitative_columns <- ordered_quantitative_columns
  
  }
  
  # Substitute zero abundances with NA.
  report$Data[ , report$Quantitative_columns][report$Data[ , report$Quantitative_columns] == 0] <- as.numeric(x = NA) 
  
  # Remove features that have missing values to all quantitative columns.
  missing_abundances_mask <- is.na(x = report$Data[ , report$Quantitative_columns])
  completely_empty_features <- apply(missing_abundances_mask, 1, all, na.rm = TRUE)
  report$Data <- report$Data[!completely_empty_features, ]
  
  # Perform scaling per fasta file.
  if(fasta_scale){
    
    major_fasta_files <- unlist(x = lapply(report$Data$fasta_file, function(x){ return(strsplit(x = x, split = ";")[[1]][1])}))
    unique_fasta_files <- unique(x = major_fasta_files)
    
    if(length(x = unique_fasta_files) > 1){
    
      per_fasta_file_ids <- lapply(unique_fasta_files, function(x){
        
        return(which(x = major_fasta_files == x))
        
      })
      
      fasta_factors <- lapply(per_fasta_file_ids, function(x){
        
        return(colSums(x = report$Data[x, report$Quantitative_columns], na.rm = TRUE))
        
      })
      
      for(i in 1:length(x = fasta_factors)){
        
        report$Data[per_fasta_file_ids[[i]], report$Quantitative_columns] <- t(x = t(x = report$Data[per_fasta_file_ids[[i]], report$Quantitative_columns])/fasta_factors[[i]])
        
      }
    
    } else {
      
      warning("Warning at function preprocessing():\n\tThere is only a single fasta file found in the dataset, per fasta file scaling won't be performed.")
    
    }
    
  }
  
  # Perform minimum quantitative value filtering per condition.
  if(quantified_per_condition > 0){
    
    existing_adundances_mask <- !is.na(x = report$Data[ , report$Quantitative_columns])
    existing_adundances_per_condition <- as.data.frame(lapply(design, function(x){
      
      return(rowSums(x = existing_adundances_mask[ , na.omit(object = x)]))
      
    }), stringsAsFactors = FALSE)

    features_to_keep <- apply(existing_adundances_per_condition >= quantified_per_condition, 1, function(x){ return(all(x, na.rm = TRUE))})
    
    report$Data <- report$Data[features_to_keep, ]
    
  }
  
  # log2 transformation.
  if(perform_log2){
    
    report$Data[ , report$Quantitative_columns] <- log2(x = report$Data[ , report$Quantitative_columns])
  
  }
    
  # Normalization.
  if(perform_normalization){

    if(normalization_method == "median"){
      
      # Calculate the median abundance per quantitative column.
      column_factors <- apply(report$Data[ , report$Quantitative_columns], 2, median, na.rm = TRUE)
      
    } else {
      
      # Calculate the average abundance per quantitative column.
      column_factors <- apply(report$Data[ , report$Quantitative_columns], 2, mean, na.rm = TRUE)
      
    }
    
    # If scaled by fasta file then zero-center normalization.
    if(fasta_scale){
      
      report$Data[ , report$Quantitative_columns] <- data.frame(t(t(report$Data[ , report$Quantitative_columns]) - column_factors), stringsAsFactors = FALSE)
      
    } else {
      
      # Global normalization factor is the average of column factors.
      global_normalization_factor <- mean(x = column_factors, na.rm = TRUE)
      
      # Perform normalization.
      report$Data[ , report$Quantitative_columns] <- data.frame(t(t(report$Data[ , report$Quantitative_columns]) - column_factors + global_normalization_factor), stringsAsFactors = FALSE)

    }
    
  }
  
  # Imputation.
  if(imputation){

    # Set a seed for reproducibility.
    set.seed(123)

    report$Data[ , report$Quantitative_columns] <- impute_downshift_gaussian(data = report$Data[ , report$Quantitative_columns],
                                                                               design = design,
                                                                               condition_sensitive = conditional_imputation)

  }
  
  # Rename rows.
  rownames(x = report$Data) <- 1:nrow(x = report$Data)
  
  return(report)
  
}

impute_downshift_gaussian <- function(data = NULL, design = NULL, condition_sensitive = FALSE, downshift = 0.4, width = 0.3, min_observations = 3) {
  
  data_imp <- data
  
  if(condition_sensitive){
  
    for(protein in 1:nrow(x = data)){
      
      for(condition in colnames(x = design)){
        
        samples <- unlist(x = na.omit(object = design[,condition]))
        protein_condition_data <- data[protein, samples]
        
        # Impute only if all replicates are missing.
        if(all(is.na(x = protein_condition_data))){

          observations <- unlist(x = data[ , samples])
          observations <- observations[!is.na(x = observations)]
          
          if(length(x = observations) >= min_observations){

            obs_mean <- mean(x = observations)
            obs_sd <- sd(x = observations)

            imputation_mean <- obs_mean - downshift * obs_sd
            imputation_sd <- width * obs_sd

            data_imp[protein, samples] <- rnorm(n = length(x = samples), mean = imputation_mean, sd = imputation_sd)
            
          }
          
        }
        
      }
      
    }
  
  } else {
    
    for(sample in colnames(x = data)){
      
      sample_data <- data[ , sample]
      sample_na <- is.na(x = sample_data)
      
      if(any(sample_na)){
        
        observations <- sample_data[!sample_na]
        
        if(length(x = observations) >= min_observations){
          
          obs_mean <- mean(x = observations)
          obs_sd <- sd(x = observations)
          
          imputation_mean <- obs_mean - downshift * obs_sd
          imputation_sd <- width * obs_sd
          
          data_imp[sample_na, sample] <- rnorm(n = sum(sample_na), mean = imputation_mean, sd = imputation_sd)
          
        }
        
      }
      
    }
    
  }

  return(data_imp)
  
}

perform_statistics <- function(report = NULL, metadata = NULL){
  
  # Create long version of the Design.
  samples <- metadata$Design %>%
    tidyr::pivot_longer(cols = tidyr::all_of(x = colnames(x = metadata$Design)), names_to = "Condition", values_to = "Sample") %>%
    tidyr::drop_na(Sample) %>%
    dplyr::mutate(Condition = factor(x = Condition, levels = colnames(x = metadata$Design))) %>%
    dplyr::arrange(Condition)

  # Create a design matrix.
  design <- model.matrix(object = ~0 + Condition, data = samples)
  colnames(x = design) <- colnames(x = metadata$Design)

  # Validate and adjust the contrast matrix.
  contrast_matrix <- metadata$Contrasts

  if(nrow(x = contrast_matrix) > 1){

    if(colnames(x = contrast_matrix)[2] == "All"){

      reference_condition <- contrast_matrix$All == 1

      if(sum(reference_condition, na.rm = TRUE) != 1){

        stop("Error at function perform_statistics():\n\tWhen All column is used in the contrast table only a single condition may be the reference (value of 1).")

      }

      if(sum(contrast_matrix$All == 2, na.rm = TRUE) != (nrow(x = contrast_matrix) - 1)){

        stop("Error at function perform_statistics():\n\tWhen All column is used in the contrast table all conditions, besides the reference one, should have the value of 2.")

      }

      # Create all the possible comparisons with the reference condition.
      comparisons <- paste0(contrast_matrix$Conditions[reference_condition], "_vs_", contrast_matrix$Conditions[!reference_condition])
      # Create a diagonal matrix denoting the conditions to be compared to the reference.
      comparison_pairs <- diag(x = rep(x = -1, times = (nrow(x = contrast_matrix) - 1)))

      # Create the contrast matrix.
      contrast_matrix <- matrix(data = NA,
                                ncol = length(x = comparisons),
                                nrow = nrow(x = contrast_matrix),
                                dimnames = list(contrast_matrix$Condition, comparisons))

      # Add the comparison memberships for the reference and comparing conditions.
      contrast_matrix[reference_condition, ] <- rep(x = 1, times = ncol(x = contrast_matrix))
      contrast_matrix[!reference_condition, ] <- comparison_pairs

    } else {

      # Name rows according to Conditions.
      rownames(x = contrast_matrix) <- contrast_matrix$Conditions

      # Remove Conditions column.
      contrast_matrix <- contrast_matrix[, -1, drop = FALSE]

      # Adjust contrasts to fit limma format. This takes account multiple condition vs multiple condition comparison too.
      contrast_matrix <- apply(contrast_matrix, 2, function(x){

        reference_condition <- x == 1
        reference_condition_count <- sum(reference_condition, na.rm = TRUE)
        comparing_condition <- x == 2
        comparing_condition_count <- sum(comparing_condition, na.rm = TRUE)

        x[reference_condition] <- 1/reference_condition_count
        x[comparing_condition] <- -1/comparing_condition_count

        return(x)

      })

      # Check for any comparisons that reference or comparing condition does not exist, and exclude them.
      comparisons_validation <- colSums(x = contrast_matrix) == 0

      if(any(!comparisons_validation, na.rm = TRUE)){

        warning(paste0("Warning at function perform_statistics():\n\tComparison(s) ",
                paste(colnames(x = contrast_matrix)[!comparisons_validation], collapse = ", "),
                " missing either reference or comparing condition, and thus will be removed."))

        contrast_matrix <- contrast_matrix[ , comparisons_validation]

      }

    }

    # Check if there are any remaining valid comparisons and do statistics.
    if(ncol(x = contrast_matrix) == 0){

      warning("Warning at function perform_statistics():\n\tNo remaining valid comparisons. No statistical analysis can be performed.")
      return(NULL)

    } else {

      limma_fit <- suppressWarnings(expr = limma::lmFit(object = report$Data[ , report$Quantitative_columns], design = design))
      contrasts_fit <- suppressWarnings(expr = limma::contrasts.fit(fit = limma_fit, contrasts = contrast_matrix))
      bayes_fit <- suppressWarnings(expr = limma::eBayes(fit = contrasts_fit))

      statistics <- lapply(colnames(x = contrast_matrix), function(x){

        stats_table <- limma::topTable(fit = bayes_fit, coef = x, number = Inf, sort.by = "none")
        stats_table <- stats_table[ , c("AveExpr", "logFC", "adj.P.Val")]
        colnames(x = stats_table) <- c("Avg", "logFC", "FDR")
        return(stats_table)

      })

      names(x = statistics) <- colnames(x = contrast_matrix)

      return(list(Design = design, Contrasts = contrast_matrix, Stats = statistics))

    }

  } else {

    warning("Warning at function perform_statistics():\n\tOnly a single condition is found in the contrast table. No statistical analysis can be performed.")
    return(NULL)

  }

}

protein_qc_plots <- function(report = NULL, metadata = NULL, perform_log2 = NULL, normalization_method = NULL, custom_condition_colors = FALSE){
  
  # Create a condition dictionary type of vector to map samples to conditions.
  condition_dict <- na.omit(object = reshape2::melt(data = metadata$Design, id.vars = NULL))
  condition_dict <- setNames(object = as.vector(x = condition_dict$variable), nm = condition_dict$value)
  
  # Should be a vector of same length as the number of conditions, with valid colors.
  validate_custom_condition_colors <- checkmate::check_flag(x = custom_condition_colors, na.ok = FALSE, null.ok = FALSE)
  use_custom_colors <- FALSE

  if(!isTRUE(x = validate_custom_condition_colors)){
    
    validate_custom_condition_colors_length <- checkmate::check_vector(x = custom_condition_colors, len = ncol(x = metadata$Design), null.ok = FALSE)
    
    if(!isTRUE(x = validate_custom_condition_colors_length)){
      
      custom_condition_colors <- FALSE
      use_custom_colors <- FALSE
      warning("Warning at function DownstreameR():\n\tInvalid custom_condition_colors value. Should have the same length as the number of conditions.")
      
    } else {
      
      use_custom_colors <- TRUE
      
    }
    
  } else {
    
    custom_condition_colors <- FALSE
    use_custom_colors <- FALSE
    
  }
  
  # Count number of missing values per sample and create a barplot.
  NAs <- data.frame(Samples = report$Quantitative_columns,
                    Count = colSums(x = is.na(x = report$Data[ , report$Quantitative_columns]), na.rm = TRUE),
                    stringsAsFactors = FALSE)
  
  NAs$Conditions <- factor(x = condition_dict[NAs$Samples], levels = colnames(x = metadata$Design))
  NAs$Samples <- factor(x = NAs$Samples, levels = NAs$Samples)
  
  plot_NAs <- ggplot2::ggplot(data = NAs, mapping = ggplot2::aes(x = Samples, y = Count, fill = Conditions)) + 
    ggplot2::geom_bar(stat = "identity", position = "dodge", ggplot2::aes(color = Samples))  +
    ggplot2::geom_text(ggplot2::aes(label = Count, y = Count), position = ggplot2::position_dodge(1), vjust = -1) +
    ggplot2::ylab(label = "Missing values count") +
    ggplot2::xlab(label = "Samples") +
    ggplot2::ggtitle(label = "Number of missing protein expressions (NA) per sample") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    ggplot2::guides(color = "none")
  
  # Distribution plot per condition.
  # Add protein_id column.
  expression_wide <- data.frame(protein_id = rownames(x = report$Data[ , report$Quantitative_columns]),
                                report$Data[ , report$Quantitative_columns],
                                stringsAsFactors = FALSE)
  
  expression_long <- reshape2::melt(expression_wide, id.vars = "protein_id")
  colnames(x = expression_long)[2:3] <- c("Samples", "Quantity")
  expression_long$Samples <- factor(x = expression_long$Samples, levels = report$Quantitative_columns)
  expression_long$Conditions <- factor(x = condition_dict[expression_long$Samples],  levels = colnames(x = metadata$Design))
  
  plot_distributions <- ggplot2::ggplot(data = na.omit(object = expression_long), mapping = ggplot2::aes(x = Quantity, color = Conditions, y = ggplot2::after_stat(scaled))) + 
    ggplot2::geom_density() +
    ggplot2::ylab(label = "Scaled density") +
    ggplot2::xlab(label = "Protein expression (log2)") +
    ggplot2::ggtitle(label = paste0("Distributions of normalized (", normalization_method,") protein expressions per condition")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  
  # Violin plot combined with box plot of protein expressions per sample.
  plot_violin <- ggplot2::ggplot(data = na.omit(object = expression_long), mapping = ggplot2::aes(x = Samples, y = Quantity)) +
    ggplot2::geom_violin(width = 1, ggplot2::aes(fill = Conditions)) +
    ggplot2::geom_boxplot(width = 0.1, ggplot2::aes(color = Samples)) +
    ggplot2::ylab(label = "Protein expression (log2)") +
    ggplot2::xlab(label = "Samples") +
    ggplot2::ggtitle(label = paste0("Violin plot of normalized (", normalization_method,") protein expressions per condition")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    ggplot2::guides(color = "none", fill = "none")
  
  # Protein feature CV per condition.
  
  # If they are on log2 scale, take abundances back to linear scale.
  if(perform_log2){
    
    expression_long$Quantity <- 2^expression_long$Quantity 
    
  }
  
  # Calculate CVs.
  cvs_long <- expression_long %>% 
              dplyr::group_by(protein_id, Conditions) %>% 
              dplyr::summarise(CV = 100*(sd(x = Quantity, na.rm = TRUE) / mean(x = Quantity, na.rm = TRUE)), .groups = "drop")
  
  # Create a CVs summary table to annotate on the plot.
  cvs_annotation <- cvs_long %>%
                    dplyr::group_by(Conditions) %>%
                    dplyr::summarise("Mean CV" = round(x = mean(CV, na.rm = TRUE), digits = 2), 
                                     "Median CV" = round(x = median(CV, na.rm = TRUE), digits = 2), .groups = "drop")
  
  plot_CV <- ggplot2::ggplot(data = na.omit(object = cvs_long), mapping = ggplot2::aes(x = Conditions, y = CV)) + 
    ggplot2::geom_boxplot(ggplot2::aes(color = Conditions)) +
    ggplot2::xlab(label = "Conditions") +
    ggplot2::ylab(label = "Protein expression CV (%)") +
    ggplot2::ggtitle(label = paste0("Per condition CV boxplot (total mean: ", round(mean(x = cvs_long$CV, na.rm = TRUE), digits = 2), "%, total median: ", round(median(x = cvs_long$CV, na.rm = TRUE), digits = 2), "%)")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    ggplot2::ylim(0, max(cvs_long$CV, na.rm = TRUE) + 20) +
    ggplot2::annotation_custom(grob = gridExtra::tableGrob(t(x = cvs_annotation[ , c(2,3)]), 
                                                           theme = gridExtra::ttheme_default(base_size = 10), 
                                                           rows = colnames(x = cvs_annotation)[2:3], 
                                                           cols = cvs_annotation$Conditions), 
                               xmin = -Inf, 
                               xmax = Inf,
                               ymin = max(cvs_long$CV, na.rm = TRUE))
  
  # Sample-wise heatmap.
  sample_correlations <- cor(x = report$Data[, report$Quantitative_columns], use = "pairwise.complete.obs", method = "pearson")
  sample_correlations <- as.data.frame(x = as.table(x = sample_correlations), stringsAsFactors = FALSE)
  colnames(x = sample_correlations) <- c("Sample_1", "Sample_2", "Correlation")
  sample_correlations$Sample_1 <- factor(x = sample_correlations$Sample_1, levels = report$Quantitative_columns)
  sample_correlations$Sample_2 <- factor(x = sample_correlations$Sample_2, levels = report$Quantitative_columns)
  
  plot_heatmap <- ggplot2::ggplot(data = sample_correlations, mapping = ggplot2::aes(x = Sample_1, y = Sample_2, fill = Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)), color = "white", size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_gradient2(low = "#1a9850", mid = "yellow", high = "red", midpoint = 0.5) +
    ggplot2::ylab(label = "Samples") +
    ggplot2::xlab(label = "Samples") +
    ggplot2::ggtitle(label = "Sample correlation heatmap") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  
  # PCA analysis and PCA plot of the two first components.
  pca <- stats::prcomp(formula = ~ .,data =  report$Data[ , report$Quantitative_columns], na.action = na.omit)
  pca_df <- data.frame(pca$rotation[, 1:2], 
                       Samples = factor(x = rownames(x = pca$rotation), levels = report$Quantitative_columns),
                       Conditions = factor(x = condition_dict[rownames(x = pca$rotation)], levels = colnames(x = metadata$Design)),
                       stringsAsFactors = F)
  
  plot_PCA <- ggplot2::ggplot(data = pca_df, mapping = ggplot2::aes(x = PC1, y = PC2, color = Conditions, fill = Samples)) +
    ggplot2::geom_point(size = 3, shape = 21, stroke = 2) +
    ggplot2::xlab(label = paste0("Component ", 1, " - ", round(x = pca$sdev[1]^2/(sum(pca$sdev^2, na.rm = TRUE)), digits = 2)*100, " %" )) +
    ggplot2::ylab(label = paste0("Component ", 2, " - ", round(x = pca$sdev[2]^2/(sum(pca$sdev^2, na.rm = TRUE)), digits = 2)*100, " %" )) + 
    ggplot2::ggtitle(label = "PCA plot of per sample relative protein expressions (Components 1 and 2)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50, 50, 50))
  
  if(use_custom_colors){
    
    plot_NAs <- plot_NAs + ggplot2::scale_fill_manual(values = setNames(object = custom_condition_colors, nm = colnames(x = metadata$Design)))
    plot_distributions <- plot_distributions + ggplot2::scale_color_manual(values = setNames(object = custom_condition_colors, nm = colnames(x = metadata$Design)))
    plot_violin <- plot_violin + ggplot2::scale_fill_manual(values = setNames(object = custom_condition_colors, nm = colnames(x = metadata$Design)))
    plot_CV <- plot_CV + ggplot2::scale_color_manual(values = setNames(object = custom_condition_colors, nm = colnames(x = metadata$Design)))
    plot_PCA <- plot_PCA + ggplot2::scale_color_manual(values = setNames(object = custom_condition_colors, nm = colnames(x = metadata$Design)))
    
  }
  
  return(list(NA_plot = plot_NAs,
              Distribution_plot = plot_distributions,
              Violin_plot = plot_violin,
              CV_plot = plot_CV,
              Heatmap = plot_heatmap,
              PCA_plot = plot_PCA))
  
}

plot_statistics <- function(statistics = NULL, fdr = 0.05, log2_fc = 0.58496250073){
  
  
  comparison_plots <- vector(mode = "list", length = length(x = statistics) + 1)
  names(x = comparison_plots) <- c("Regulation_counts", names(x = statistics))
  
  page_count <- 4
  
  for(comparison in names(x = statistics)){
    
    comparison_df <- statistics[[comparison]]
    comparison_df$Regulation <- factor(x = comparison_df$Regulation, levels = c("Up-regulated", "Down-regulated", "Non-regulated"))
    
    plot_volcano <- ggplot2::ggplot(data = na.omit(object = comparison_df), mapping = ggplot2::aes(x = logFC, y = -log10(x = FDR), color = Regulation)) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::geom_vline(xintercept = c(-log2_fc, log2_fc), linetype = "dashed", color = "black") +
      ggplot2::geom_hline(yintercept = -log10(x = fdr), linetype = "dashed", color = "black") +
      ggplot2::annotate(geom = "text", x = 1.05, y = max(comparison_df$FDR, na.rm = TRUE) * 0.95, label = paste0("log2FC = ", round(x = log2_fc, digits = 3)), hjust = 0, size = 3) +
      ggplot2::annotate(geom = "text", x = -1.05, y = max(comparison_df$FDR, na.rm = TRUE) * 0.95, label = paste0("log2FC = ", round(x = -log2_fc, digits = 3)), hjust = 1, size = 3) +
      ggplot2::annotate(geom = "text", x = max(comparison_df$logFC, na.rm =  TRUE) * 0.8, y = -log10(0.05) + 0.1, label = "FDR = 0.05", vjust = 0, size = 3) +
      ggplot2::scale_color_manual(values = c("Up-regulated" = "#80ff80", "Down-regulated" = "#ff4d4d", "Non-regulated" = "grey")) +
      ggplot2::xlab(label = "log2 Fold Change") +
      ggplot2::ylab(label = "-log10 FDR") +
      ggplot2::ggtitle(label = paste0("Volcano plot of ", comparison)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50))
    
    plot_MA <- ggplot2::ggplot(data = na.omit(object = comparison_df), mapping = ggplot2::aes(x = Avg, y = logFC, color = Regulation)) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::scale_color_manual(values = c("Up-regulated" = "#80ff80", "Down-regulated" = "#ff4d4d", "Non-regulated" = "grey")) +
      ggplot2::xlab(label = "Log2 Average Expression") +
      ggplot2::ylab(label = "Log2 Fold Change") +
      ggplot2::ggtitle(label = paste0("MA plot of ", comparison)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50))
    
    comparison_plots[[comparison]] <- (plot_volcano / plot_MA) + 
      patchwork::plot_annotation(title = "Differential Expression Analysis",
                                 subtitle = paste0("Comparison ", comparison), 
                                 caption = paste0("Page ", page_count),
                                 theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50)))
    
    page_count <- page_count + 1
    
  }
  
  # Do the protein regulation barplot.
  regulation_counts <- as.data.frame(x = do.call("rbind", lapply(statistics, function(x){
    
    return(table(factor(x = x$Regulation, levels = c("Up-regulated", "Down-regulated", "Non-regulated"))))
    
  })), stringsAsFactors = FALSE)
    
  regulation_counts$Comparison <- names(x = statistics)
  regulation_counts <- reshape2::melt(data = regulation_counts, id.vars = "Comparison")
  regulation_counts$Comparison <- factor(x = regulation_counts$Comparison, levels = names(x = statistics))
  colnames(x = regulation_counts)[2:3] <- c("Regulation", "Count")
  
  plot_count_bars <- ggplot2::ggplot(data = regulation_counts, mapping = ggplot2::aes(x = Comparison, y = Count, fill = Regulation)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5) +
    ggplot2::geom_text(data = subset(x = regulation_counts, subset = Count != 0), mapping = ggplot2::aes(label = Count), position = ggplot2::position_stack(vjust = 0.7), color = "white", size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("Up-regulated" = "#80ff80", "Down-regulated" = "#ff4d4d", "Non-regulated" = "grey")) +
    ggplot2::ylab("Protein Count") +
    ggplot2::xlab("Comparison") +
    ggplot2::ggtitle("Protein regulation barplot") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  
  comparison_plots[["Regulation_counts"]] <- (plot_count_bars / (ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50)))) + 
    patchwork::plot_annotation(title = "Differential Expression Analysis", subtitle = "Overview", theme = ggplot2::theme(plot.margin = ggplot2::margin(t = 50, b = 50))) + 
    patchwork::plot_layout(heights = c(1, 1))
  
  return(comparison_plots)
  
}
