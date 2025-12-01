# Set base directory
base_dir <- "/data/users/david/hungate_410_asmash"

# Load required libraries
library(ggtree)
library(scales)
library(ggsci)
library(ggtreeExtra)
library(treeio)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(rcartocolor)
library(forcats)

# Read phylogenetic tree
tree <- read.tree(file.path(base_dir, "../hungate_genomes/assign_tax_to_hun/phylophlan/hun410/hun410_rooted.treefile"))

# Load taxonomy data
tax_file <- fread(file.path(base_dir, "../hungate_genomes/assign_tax_to_hun/ani/gtdbtk.ani_closest.tsv"), fill = TRUE)
tax <- data.frame(tax_file, do.call(rbind, str_split(tax_file$reference_taxonomy, ";")))
colnames(tax)[7:13] <- c("domain", "phylum", "class", "order", "family", "genus", "species")
colnames(tax)[1] <- "genome_id"
setDT(tax)
tax[genome_id == "MetmilDSM16643", phylum := "p__Methanobacteriota"]
tax[genome_id == "ButspAE2032", phylum := "p__Firmicutes_A"]

# Load product data
asmash_product_list <- fread(file.path(base_dir, "antismash_hungate_410_product_list.tsv")) %>%
  mutate(contig = sub("\\.region.*.", "", filename)) %>%
  rename(genome_id = folder_name, bgc_type = product) %>%
  select(genome_id, bgc_type)

ribosomal_cols <- c(
  "RaS-RiPP", "RRE-containing", "RiPP-like", "lanthipeptide-class-i", "lanthipeptide-class-ii",
  "cyclic-lactone-autoinducer", "lanthipeptide-class-iv", "lassopeptide", "thiopeptide", "LAP",
  "sactipeptide", "ranthipeptide", "proteusin"
)

product_tree <- asmash_product_list %>%
  filter(bgc_type %in% ribosomal_cols)

# Transform data to wide format
wide_df <- product_tree %>%
  distinct() %>%
  mutate(value = bgc_type) %>%
  pivot_wider(
    id_cols = genome_id,
    names_from = bgc_type,
    values_from = value,
    values_fill = list(value = NA)
  )

annot_df <- column_to_rownames(wide_df, var = "genome_id")

# Taxonomy annotation
annot_tax_df <- tax %>%
  select(genome_id, phylum) %>%
  rename(label = genome_id) %>%
  distinct()

# Create phylogenetic tree with annotations
circ <- ggtree(tree, layout = "fan", open.angle = 90)
circ2 <- circ %<+% annot_tax_df +
  geom_tippoint(pch = 16, aes(col = phylum), size = 1) +
  scale_color_manual(values = carto_pal("ag_GrnYl", n = n_distinct(annot_tax_df$phylum)))

# Create heatmap
p1 <- gheatmap(
  circ2, annot_df, legend_title = "BGC type", custom_column_labels = "",
  offset = 0, width = 0.4, color = "#cccccc98",
  colnames_angle = 0, colnames_offset_y = 0,
  font.size = 2, colnames_offset_x = 0
) +
  scale_fill_manual(
    values = c(pal_jco("default")(10), pal_jco("default", alpha = 0.4)(4)),
    na.value = "white"
  )

p2 <- rotate_tree(p1, +90) +
  theme_tree(legend.position = c(0, 0)) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -50, unit = "pt"),
    legend.text = element_text(size = 6, face = "bold", family = "sans"),
    legend.title = element_text(size = 6, face = "bold", family = "sans")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

# Save tree plots
ggsave(p2, filename = file.path(base_dir, "analysis/figures/phylogenetic_tree_with_ribosomal_bgc.png"), width = 7, height = 7, units = "in", dpi = 300)
ggsave(p2, filename = file.path(base_dir, "analysis/figures/phylogenetic_tree_with_ribosomal_bgc.pdf"), width = 7, height = 7, units = "in", dpi = 300)
ggsave(p2, filename = file.path(base_dir, "analysis/figures/phylogenetic_tree_with_ribosomal_bgc.svg"), width = 7, height = 7, units = "in", dpi = 300)

# Define colors for BGC types
bgc_colors <- c(
  "ranthipeptide" = "#3B3B3BFF",
  "lanthipeptide-class-ii" = "#868686FF",
  "cyclic-lactone-autoinducer" = "#0073C2FF",
  "RRE-containing" = "#7AA6DCFF",
  "RiPP-like" = "#4A6990FF",
  "lassopeptide" = "#003C67FF",
  "thiopeptide" = "#86868666",
  "RaS-RiPP" = "#A73030FF",
  "lanthipeptide-class-i" = "#EFC000FF",
  "LAP" = "#7AA6DCFF",
  "proteusin" = "#8F7700FF",
  "sactipeptide" = "#EFC00066",
  "lanthipeptide-class-iv" = "#CD534CFF"
)

# Create bar plot for product counts
product_count_p1 <- product_tree %>%
  ggplot(aes(x = fct_rev(fct_infreq(bgc_type)))) +
  geom_bar(aes(fill = bgc_type)) +
  coord_flip() +
  theme_classic() +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = bgc_colors) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8, face = "bold", family = "sans"),
    axis.text.x = element_text(size = 5, face = "bold", family = "sans")
  )

# Save bar plot
ggsave(product_count_p1, filename = file.path(base_dir, "analysis/figures/ribosomal_bgc_count_bar.png"), width = 3, height = 3, units = "in", dpi = 300)
ggsave(product_count_p1, filename = file.path(base_dir, "analysis/figures/ribosomal_bgc_count_bar.pdf"), width = 5, height = 5, units = "in", dpi = 300)
