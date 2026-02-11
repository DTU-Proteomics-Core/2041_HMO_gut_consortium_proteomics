setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
library(patchwork)

stoichiometries <- read.csv(file = "Data/iBAQ_stoichiometries.csv", header = TRUE, sep = ",")

# Add the Roseburia species summation.
roseburia_sum <- stoichiometries %>%
                  dplyr::filter(grepl(pattern = "^Roseburia", x = Species)) %>%
                  dplyr::group_by(Sample) %>%
                  dplyr::summarise(Count = sum(Count), .groups = "drop") %>%
                  dplyr::mutate(Species = "Roseburia species") %>%
                  dplyr::select(Species, Sample, Count)

stoichiometries <- rbind(stoichiometries, roseburia_sum)

# Add condition levels.
stoichiometries$Condition <- sapply(as.character(x = stoichiometries$Sample), function(x){ return(strsplit(x = x, split = "_")[[1]][2])})
stoichiometries$Condition <- factor(x = stoichiometries$Condition, levels = c("HMOs", "Mix", "Fiber"))

# Define comparisons.
comparisons <- list(c("Mix", "HMOs"), c("Mix", "Fiber"), c("HMOs", "Fiber"))

# Color palette
colors <- c("#C03830", "#317EC2", "#008835")

# Function to create individual panels.
plot_panel <- function(data = NULL, comparisons = NULL, palette = NULL){

  # Custom function to calculate and render p-value labels in scientific format.
  add_pval_labels <- function(data, comparisons) {
    
    format_scientific <- function(p) {
      
      if(p < 0.01){
        
        exponent_value <- floor(x = log10(x = p))
        mantissa <- p / 10^exponent_value
        exponent_abs <- abs(x = exponent_value)
        superscript_digits <- c("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
        exponent_char <- strsplit(x = as.character(x = exponent_abs), "")[[1]]
        exponent_superscript <- paste(superscript_digits[as.numeric(x = exponent_char) + 1], collapse = "")
        
        return(paste0("P = ", sprintf("%.1f", mantissa), " × 10⁻", exponent_superscript))
        
      } else {
        
        return(paste0("P = ", sprintf("%.3f", p)))
      
      }
    }
    
    stat_test <- data %>% rstatix::t_test(formula = Count ~ Condition, comparisons = comparisons) %>%
                          rstatix::add_xy_position(x = "Condition")
    
    stat_test$p_formatted <- sapply(stat_test$p, function(x){return(format_scientific(p = x))})
    
    return(stat_test)
    
  }
  
  # Calculate stats and form custom labels.
  stat_test <- add_pval_labels(data = data,comparisons =  comparisons)
  
  figure <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = Condition, y = Count, fill = Condition)) +
              ggplot2::stat_summary(fun = mean, geom = "col", width = 0.5) +
              ggplot2::stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
              ggplot2::geom_point(size = 2, shape = 16) + 
              ggpubr::stat_pvalue_manual(data = stat_test, label = "p_formatted", size = 3) +
              ggplot2::facet_wrap(facets = ~ Species, nrow = 1) +
              ggplot2::scale_fill_manual(values = palette) +
              ggplot2::theme_classic() +
              ggplot2::theme(legend.position = "none", 
                             strip.background = element_blank(),
                             strip.text = element_text(size = 10)) +
              ggplot2::labs(y = "Relative abundance (%)", x = "Condition")


  return(figure)

}

panel_species <- unique(stoichiometries$Species)[c(1:2, 6, 3:5)]
panel_plots <- lapply(panel_species, function(x){
  
  return(plot_panel(data = stoichiometries[stoichiometries$Species == x, ], comparisons = comparisons, palette = colors))
  
})

combined_plot <- patchwork::wrap_plots(panel_plots, ncol = 3, nrow = 2) + patchwork::plot_annotation(tag_levels = 'A')


ggplot2::ggsave("Data/Supplementary_figure.pdf", plot = combined_plot, width = 15, height = 10, units = "in", device = cairo_pdf)

