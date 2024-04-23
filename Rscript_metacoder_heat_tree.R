#Metacoder


library(metacoder)
library(cowplot)
library(ggplot2)
library(microViz)

setwd = ("C/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_16S/All16S/Metacoder/Figure_metacoder_article") 


hmp_otus <- read.csv2("C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_16S/All16S/Metacoder/Figure_metacoder_article/Metacoder_18S_PSRosup10.csv", sep=";", header=TRUE)
hmp_otus

hmp_samples <- read.csv2("C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_16S/All16S/Metacoder/Figure_metacoder_article/PSRo_samples.csv", sep=";", header=TRUE)
hmp_samples


tm_obj = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
                        class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                        class_regex = "^(.+)__(.+)$")


# Convert counts to proportions
tm_obj$data$otu_table <- calc_obs_props(tm_obj, data = "tax_data", cols = hmp_samples$sample_id)

# Get per-taxon counts
tm_obj$data$tax_table <- calc_taxon_abund(tm_obj, data = "otu_table", cols = hmp_samples$sample_id)

tm_obj

hmp_samples
unique(hmp_samples$Health)

comparisons <- list(c("Normal", "Disease"))

# Calculate difference between groups
tm_obj$data$diff_table <- compare_groups(tm_obj, data = "tax_table",
                                         cols = hmp_samples$sample_id,
                                         groups = hmp_samples$Health,
                                         combinations = comparisons)
tm_obj$data$diff_table

## Make comparison plots

plot_comp <- function(comp_pair) {
  set.seed(1) # so plot layout and sizes always look the same
  tm_obj %>%
    filter_obs(data = "diff_table", treatment_1 %in% comp_pair & treatment_2 %in% comp_pair) %>%
    heat_tree(node_size = n_obs,
              node_size_range = c(0.01, 0.05),
              node_color = log2_median_ratio,
              node_label = taxon_names,
              node_color_range = diverging_palette(),
              node_color_trans = "linear",
              node_color_interval = c(-3, 3),
              edge_color_interval = c(-3, 3),
              node_size_axis_label = "Number of ASV",
              node_color_axis_label = "Log2 ratio median proportions",
              title = paste0(comp_pair[1], ' vs. ', comp_pair[2]),
              make_node_legend = TRUE, 
              make_edge_legend = TRUE)
}


comp_plots <- lapply(comparisons, plot_comp)
comp_plots

#Colors red for diseases blue for healthy, visualiziation at species level, but remove NA names in species
#heat tree with different layout; up to species level

tm_obj %>% 
  filter_taxa(taxon_ranks == "s", supertaxa = TRUE) %>% # subset to the species rank
  filter_taxa(taxon_ranks != "s" | !grepl("NA", taxon_names)) %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_size_range = c(0.01, 0.05),
            node_color = log2_median_ratio,
            node_color_axis_label = "Log2 ratio median proportions",
            node_color_interval = c(-3, 3),
            edge_color_interval = c(-3, 3),
            node_color_range = (c("#D88C9A",'#E5DCC5',"#81C4C0")),
            layout = "davidson-harel", initial_layout = "reingold-tilford")

#heat tree with different layout; up to genus level (removing NA)

tm_obj %>% 
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the species rank
  filter_taxa(taxon_ranks != "g" | !grepl("NA", taxon_names)) %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_size_range = c(0.01, 0.05),
            node_color = log2_median_ratio,
            node_color_axis_label = "Log2 ratio median proportions",
            node_color_interval = c(-3, 3),
            edge_color_interval = c(-3, 3),
            node_color_range = (c("#D88C9A",'#E5DCC5',"#81C4C0")),
            layout = "davidson-harel", initial_layout = "reingold-tilford")
