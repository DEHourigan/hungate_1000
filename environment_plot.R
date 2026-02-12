########################################
## Packages
########################################
.libPaths(c("/data/san/data0/users/david/rstudio/packages", .libPaths()))
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(ape)
library(phangorn)
library(ggtree)
library(stringr)
library(tidytable)
library(janitor)
library(ggtreeExtra)
library(purrr)
library(forcats)
library(thacklr)
library(gggenomes)
library(patchwork)

set.seed(1738)
setwd("/data/san/data3/users/david/hungate")

########################################
habitat_df <- fread("data/habitat.out.smorfs.tsv")
protein_acc_df <- fread("data/corepeptide_diamond.nr.out") %>%
  group_by(V1) %>%
  filter(V3 == 100) %>%	
  slice_min(order_by = V11, n = 1) %>%
  ungroup() %>%
  select(V1, V2) %>%
  dplyr::rename(qseqid = V1, protein_acc = V2)

number_of_core <- 513 


# plot the distribution of habitats
habitat_long <- habitat_df %>%
  separate_rows(habitat, sep = ",") %>%
  mutate(habitat = str_trim(habitat)) %>%
  left_join(protein_acc_df, by = c("qseqid"))

habitat_plot <- habitat_long %>%
	count(habitat, sort = TRUE) %>%
	mutate(is_rumen = str_detect(habitat, regex("rumen|cattle gut|goat gut|deer gut|goat rumen", ignore_case = TRUE))) %>%
	ggplot(aes(x = reorder(habitat, n), y = n, fill = is_rumen)) +
	geom_col(color = "black") +
	scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = rgb(154, 47, 45, maxColorValue = 255)), guide = "none") +
	coord_flip() +
	labs(x = "Habitat", y = "Count") +
	theme_bw()

ggsave("reports/figures/habitat_distribution.png",habitat_plot,  
	width = 15, height = 15, units = "cm", dpi = 300)

sheets <- readxl::excel_sheets("data/Supplementary_table_4.xlsx")

product_df <- purrr::map_dfr(sheets, function(sh) {
	readxl::read_excel("data/Supplementary_table_4.xlsx", sheet = sh, col_names = FALSE) %>%
		rename(col = "...1") %>%
		filter(grepl("^>", col)) %>%
		mutate(col = sub("^>", "", col)) %>%
		tidyr::separate(col, into = c("qseqid", "product"), sep = " ", extra = "merge", fill = "right") %>%
		mutate(bgc_class = sh) %>%
		select(bgc_class, qseqid, product)
}) 


########################################
## summary stats on classes
########################################
region_df <- readxl::read_excel("data/Supplementary_table_4.xlsx", sheet = 1, col_names = TRUE) %>%
	clean_names() 


region_df %>%
	group_by(bgc_type) %>%
	summarise(
		n = n(),
		n_with_region = sum(!is.na(region_id)),
		total = n_distinct(region_df$region_id)
	)

region_df %>%
	group_by(bgc_type) %>%
	filter(newly_discovered == "Yes") %>%
	summarise(
		n = n(),
		n_with_region = sum(!is.na(region_id)),
		total = n_distinct(region_df$region_id)
	)
	
# presence / absence of protein_acc by habitat, faceted by bgc_class
pa_df <- habitat_product %>%
	distinct(bgc_class, habitat, protein_acc) %>%
	mutate(present = 1L) %>%
	group_by(bgc_class) %>%
	tidyr::complete(
		habitat = unique(habitat),
		protein_acc = unique(protein_acc),
		fill = list(present = 0L)
	) %>%
	ungroup() %>%
	filter(bgc_class %in% c("Class II lanthipeptides", "Class I lanthipeptides",
		"Lassopeptides", "Ranthipeptide", "Core Peptides IIA", "Core Peptides IIC")) %>%
	filter(!is.na(protein_acc))

fwrite(pa_df, "reports/bacteriocin_core_presence_absence.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

heat_df <- habitat_product %>%
	distinct(bgc_class, habitat, protein_acc) %>%
	count(bgc_class, habitat, name = "n")

# total unique protein_acc per bgc_class
totals <- habitat_product %>%
	distinct(bgc_class, protein_acc) %>%
	count(bgc_class, name = "total_in_class")

	# compute proportion = n / total_in_class
heat_df_norm <- heat_df %>%
	left_join(totals, by = "bgc_class") %>%
	mutate(prop = ifelse(total_in_class > 0, n / total_in_class, 0)) %>%
	ungroup()

gut_rumen_df <- heat_df_norm %>%
	filter(grepl("gut|rumen", habitat, ignore.case = TRUE)) %>%
	filter(total_in_class >= 3)

heatmap <- ggplot(gut_rumen_df, aes(x = habitat, y = bgc_class, fill = prop)) +
	geom_tile(color = "grey80") +
	labs(x = "BGC Class", y = "Habitat", fill = "Proportion\n(of class)") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	# highlight "rumen"
	geom_tile(data = subset(gut_rumen_df, 
		grepl("rumen|cattle gut|goat gut|deer gut", 
		habitat, ignore.case = TRUE)),
		color = "black", size = 1, fill = NA) +
# scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) 
	scale_fill_gradient(
		low = rgb(231, 185, 42, maxColorValue = 255),
		high = rgb(154, 47, 45, maxColorValue = 255),
		labels = scales::percent_format(accuracy = 1)
	)


# red  154	47	45	
# yellow  231	185	42	
# green

ggsave(heatmap, 
  filename = "reports/figures/heatmap_bgc_habitat.png", 
  width = 20, height = 7, dpi = 300, units="cm")

fwrite(gut_rumen_df, "reports/proportions_bgc_habitat_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)





########################################
# read in alignment.out.smorfs.tsv in data/ 
########################################

alignments <- fread("data/alignment.out.smorfs.tsv") %>% 
	left_join(protein_acc_df, by = c("V1" = "qseqid")) %>%
	select(protein_acc, V1, V2, V3,V4,V13) %>%
	dplyr::rename(
		qseqid = V1, 
		target = V2, 
		seq_q = V3, 
		seq_t = V4, 
		pident = V13) %>%
	distinct()

gmsc_samples <- fread("data/100AA_gmsc_hits_hungate_core_with_samples.tsv")

# make sample long format (default the field is comma separated)
gmsc_samples_long <- gmsc_samples %>%
	separate_rows(sample, sep = ",") %>%
	mutate(sample_accession = str_trim(sample))

meta_gmsc <- fread("/data/san/data3/databases/gmsc/100/GMSC10.metadata.tsv")

product_df
protein_acc_df

gmsc_samples_long_metadata <- left_join(gmsc_samples_long, meta_gmsc, 
	by = "sample_accession") %>%
	left_join(., alignments, by = c("target")) %>%
	left_join(., product_df, by = c("qseqid"))



########################################
# 100% ID only
########################################

gmsc_samples_long_metadata_100 <- gmsc_samples_long_metadata %>%
	filter(pident == 100)

fwrite(gmsc_samples_long_metadata_100, 
	"reports/gmsc_samples_long_metadata_100ID.tsv", sep = "\t")

library(reactable)

colnames(gmsc_samples_long_metadata_100)


gmsc_samples_long_metadata_100 %>%
	filter(grepl("class II lanthipeptide", bgc_class, ignore.case = TRUE)) %>%
	group_by(qseqid, general_envo_name) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) %>% View


gmsc_samples_long_metadata_100 %>%
	filter(grepl("class II lanthipeptide", bgc_class, ignore.case = TRUE)) %>%
	filter(grepl("LJCFBI", qseqid, ignore.case = TRUE)) %>%
	group_by(general_envo_name) %>%
	select(general_envo_name,sample_accession, geographic_location) %>%
	distinct() %>%
	summarise(n = n()) 

gmsc_samples_long_metadata_100 %>%
	filter(grepl("class II lanthipeptide", bgc_class, ignore.case = TRUE)) %>%
	filter(grepl("cattle",general_envo_name, ignore.case = TRUE)) %>%
	group_by(general_envo_name) %>%
	select(general_envo_name,sample_accession, geographic_location,qseqid) %>%
	distinct() %>% 
	filter(grepl("DDFDAA", qseqid, ignore.case = TRUE)) 

unique(gmsc_samples_long_metadata_100$bgc_class)

make_faceted_plot <- function(class_name, filename, log_scale = FALSE, width = 14, height = 12) {
	class_df <- gmsc_samples_long_metadata_100 %>%
		filter(bgc_class == class_name)

	if (nrow(class_df) == 0) return(invisible(NULL))

	plot_data <- bind_rows(
		class_df %>%
			mutate(
				category = replace_na(general_envo_name, "unknown"),
				facet = "Habitat",
				fill_group = if_else(
					str_detect(category, regex("rumen|cattle gut|goat gut|deer gut|goat rumen", ignore_case = TRUE)),
					"rumen",
					"other"
				)
			) %>%
			count(facet, category, fill_group, name = "n"),
		class_df %>%
			mutate(
				category = replace_na(geographic_location, "unknown"),
				facet = "Country",
				fill_group = "country"
			) %>%
			count(facet, category, fill_group, name = "n")
	) %>%
		group_by(facet) %>%
		mutate(category = fct_reorder(category, n)) %>%
		ungroup()

	plot <- ggplot(plot_data, aes(x = category, y = n, fill = fill_group)) +
		geom_col(color = "black") +
		coord_flip() +
		facet_wrap(~ facet, scales = "free_y") +
		scale_fill_manual(
			values = c(
				"other" = "grey80",
				"rumen" = rgb(154, 47, 45, maxColorValue = 255),
				"country" = rgb(154, 47, 45, maxColorValue = 255)
			),
			guide = "none"
		) +
		labs(x = "Category", y = "Count") +
		theme_bw()

	if (log_scale) {
		plot <- plot + scale_y_log10()
	}

	ggsave(
		filename = filename,
		plot = plot,
		width = width,
		height = height,
		units = "cm",
		dpi = 300
	)
}

plot_specs <- tibble::tibble(
	bgc_class = c(
		"Class I lanthipeptides",
		"Class II lanthipeptides",
		"Core Peptides IIC",
		"Core Peptides IIA",
		"Lassopeptides"
	),
	file = c(
		"reports/figures/lanthipeptide_class_I_habitat_country_distribution.png",
		"reports/figures/lanthipeptide_class_II_habitat_country_distribution.png",
		"reports/figures/corepeptide_IIC_habitat_country_distribution.png",
		"reports/figures/corepeptide_IIA_habitat_country_distribution.png",
		"reports/figures/lassopeptide_habitat_country_distribution.png"
	),
	log_scale = c(FALSE, TRUE, FALSE, TRUE, TRUE),
	width = c(14, 18, 14, 14, 28),
	height = c(15, 15, 6, 12, 12)
)

plot_specs %>%
	pwalk(~ make_faceted_plot(..1, ..2, ..3, ..4, ..5))




########################################
## summary on how many samples each class is in, and how many unique qseqid are in each class
########################################

gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class) %>%
	summarise(
		n_samples = n_distinct(sample_accession),
		n_qseqid = n_distinct(qseqid)
	) %>%
	arrange(desc(n_samples))

# how many countries each class is in?
gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class) %>%
	summarise(n_countries = n_distinct(geographic_location)) %>%
	arrange(desc(n_countries))

# how many biomes each class is in?
gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class) %>%
	summarise(n_biomes = n_distinct(general_envo_name)) %>%
	arrange(desc(n_biomes))


########################################
# manually replot class II lanthipeptides 
########################################
class_II_df <- gmsc_samples_long_metadata_100 %>%
	filter(bgc_class == "Class II lanthipeptides")

gmsc_samples_long_metadata_100 %>%
	filter(bgc_class == "Class II lanthipeptides") %>%
	group_by(general_envo_name) %>%
	summarise(n = n()) %>%
	arrange(desc(n))  %>%
	print(n=20)


class_II_lan_stats <- gmsc_samples_long_metadata_100 %>%
	filter(bgc_class == "Class II lanthipeptides") %>%
	group_by(qseqid, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))
# which qseqid was in the most samples?
gmsc_samples_long_metadata_100 %>%
	filter(bgc_class == "Class II lanthipeptides") %>%
	group_by(qseqid) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) %>%
	head(20)

# which qseqid was in the most samples?
gmsc_samples_long_metadata_100 %>%
	filter(bgc_class == "Class II lanthipeptides") %>%
	group_by(sample_accession) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) 


# which qseqid was in the most samples?
class_II_lan_stats %>%
	filter(grepl("GEBGMI_07495", qseqid)) 

class_II_plot_data <- bind_rows(
	class_II_df %>%
		mutate(
			category = replace_na(general_envo_name, "unknown"),
			facet = "Habitat",
			fill_group = if_else(
				str_detect(category, regex("rumen|cattle gut|goat gut|deer gut|goat rumen", ignore_case = TRUE)),
				"rumen",
				"other"
			)
		) %>%
		count(facet, category, fill_group, name = "n"),
	class_II_df %>%
		mutate(
			category = replace_na(geographic_location, "unknown"),
			facet = "Country",
			fill_group = "country"
		) %>%
		count(facet, category, fill_group, name = "n")
) %>%
	group_by(facet) %>%
	mutate(category = fct_reorder(category, n)) %>%
	ungroup() %>%
	filter(n > 1)

class_II_plot <- ggplot(class_II_plot_data, aes(x = category, y = n, fill = fill_group)) +
	geom_col(color = "black") +
	coord_flip() +
	facet_wrap(~ facet, scales = "free_y", ncol = 2) +
	scale_fill_manual(
		values = c(
			"other" = "grey80",
			"rumen" = rgb(154, 47, 45, maxColorValue = 255),
			"country" = rgb(154, 47, 45, maxColorValue = 255)
		),
		guide = "none"
	) +
	labs(x = "Category", y = "Count") +
	theme_bw() +
	scale_y_log10()

ggsave(
	filename = "reports/figures/lanthipeptide_class_II_habitat_country_distribution_log.png",
	plot = class_II_plot,
	width = 28,
	height = 14,
	units = "cm",
	dpi = 300
)

########################################
# Pneumolancidin stats
########################################

pneumolancidin_country <- class_II_df %>%
	filter(grepl("LJCFBI", qseqid)) %>%
	group_by(qseqid,protein_acc, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

pneumolancidin_biome <- class_II_df %>%
	filter(grepl("LJCFBI", qseqid)) %>%
	group_by(qseqid,protein_acc, general_envo_name) %>%
	summarise(n = n()) %>%
	arrange(desc(n))





########################################
# nisin only
########################################
nisin_df <- gmsc_samples_long_metadata_100 %>%
	filter(grepl("CTPGC", seq_q, ignore.case = TRUE))

nisin_df <- bind_rows(
	nisin_df %>%
		mutate(
			category = replace_na(general_envo_name, "unknown"),
			facet = "Habitat",
			fill_group = if_else(
				str_detect(category, regex("rumen|cattle gut|goat gut|deer gut|goat rumen", ignore_case = TRUE)),
				"rumen",
				"other"
			)
		) %>%
		count(facet, category, fill_group, name = "n"),
	nisin_df %>%
		mutate(
			category = replace_na(geographic_location, "unknown"),
			facet = "Country",
			fill_group = "country"
		) %>%
		count(facet, category, fill_group, name = "n")
) %>%
	group_by(facet) %>%
	mutate(category = fct_reorder(category, n)) %>%
	ungroup() 


nisin_plot <- ggplot(nisin_df, 	aes(x = category, y = n, fill = fill_group)) +
	geom_col(color = "black") +
	coord_flip() +
	facet_wrap(~ facet, scales = "free_y", ncol = 1) +
	scale_fill_manual(
		values = c(
			"other" = "grey80",
			"rumen" = rgb(154, 47, 45, maxColorValue = 255),
			"country" = rgb(154, 47, 45, maxColorValue = 255)
		),
		guide = "none"
	) +
	labs(x = "Category", y = "Count") +
	theme_bw() +
	scale_y_log10()

ggsave(
	"reports/figures/nisin_habitat_distribution.png",
	nisin_plot,
	width = 8,
	height = 16,
	units = "cm",
	dpi = 300
)


########################################
# Pseudobutryvibrio nisin
########################################

pseudobutryvibrio_nisin <- gmsc_samples_long_metadata_100 %>%
	filter(grepl("class I lanthipeptide", bgc_class, ignore.case = TRUE)) %>%
	filter(grepl("Pseudobutyrivibrio", product)) %>%
	group_by(qseqid, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) 

n_distinct(pseudobutryvibrio_nisin$geographic_location)


########################################
# Streptococcus nisin
########################################

streptococcus_nisin <- gmsc_samples_long_metadata_100 %>%
	filter(grepl("class I lanthipeptide", bgc_class, ignore.case = TRUE)) %>%
	filter(grepl("Streptococcus", product)) %>%
	group_by(qseqid, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

########################################
# Peptostreptococcus nisin
########################################

peptostreptococcus_nisin <- gmsc_samples_long_metadata_100 %>%
	filter(grepl("class I lanthipeptide", bgc_class, ignore.case = TRUE)) %>%
	filter(grepl("JJAFAM", product)) %>%
	group_by(qseqid, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

########################################
# pediocin stats	 IIa
########################################
pediocin_set <- c(
	"LHAIJF_07965",
	"FEJLPI_14705",
	"NLCFLM_27610",
	"GILPNI_10090"
)

gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% pediocin_set) %>%
	group_by(qseqid, 
		general_envo_name, 
		geographic_location
		) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

########################################
# class IIB
########################################

class_IIB_set <- c(
  "CCKMEK_09845","CCKMEK_09850","CECIIA_09960","CECIIA_09965","MFGIFH_07670",
  "MFGIFH_07680","MFGIFH_07685","ILLMLJ_06645","ILLMLJ_06650","LHAIJF_07945",
  "LHAIJF_07950","ACJJKH_05300","ACJJKH_05305","APNHLP_17160","APNHLP_17165",
  "DDIDDG_02250","DDIDDG_02255","EGPFLP_00955","EGPFLP_00960","IKMHFE_06975",
  "IKMHFE_06980","OLMONE_04595","OLMONE_04600","ECDHJF_06265","ECDHJF_06270",
  "ECDHJF_06275","ECDHJF_06280","NOCGMI_06625","NOCGMI_06630","BKANDP_05990",
  "BKANDP_05995","BKANDP_06000","BKANDP_06005","HDEPEA_06120","HDEPEA_06125",
  "NDNJMB_02865","NDNJMB_02870","NDNJMB_09565","NDNJMB_09570","GINLMA_03675",
  "GINLMA_03680","GINLMA_04435","GINLMA_04440","GINLMA_04445","GINLMA_04450",
  "JHGBMJ_09340","JHGBMJ_09345","JHGBMJ_09350","MOCIPG_03375","MOCIPG_03380",
  "KALHLD_05605","KALHLD_05610","KCLKHH_07005","KCLKHH_07010"
)


gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% class_IIB_set) %>%
	group_by(qseqid, 
		general_envo_name, 
		# geographic_location
		) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% class_IIB_set) %>%
	group_by(qseqid, 
		# general_envo_name, 
		geographic_location
		) %>%
	summarise(n = n()) %>%
	arrange(desc(n))


gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% c("ILLMLJ_06650", "ILLMLJ_06645")) %>%
	group_by(qseqid, general_envo_name) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% c("ILLMLJ_06650", "ILLMLJ_06645")) %>%
	group_by(qseqid, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

########################################
# Ranthipeptide stats
########################################

gmsc_samples_long_metadata_100 %>%
	filter(bgc_class %in% c("Ranthipeptide")) %>%
	group_by(qseqid) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

########################################
# lassopeptide
########################################

lassopeptide_set <- unique(c(
	"GKGCKK_07220",
	"EGPPOI_06395",
	"HMHFPH_06130",
	"PHFKJA_13800",
	"GLGEAB_00345",
	"CKFMIE_00455",
	"CKFMIE_00465",
	"CKFMIE_00470",
	"MOBKGA_10540",
	"BAPACM_11765",
	"PPECAL_02995",
	"GMEKKH_03485",
	"LBGLEE_12210",
	"POBEHB_14915",
	"DLHIHL_13420"
))

lasso_biome <- gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% lassopeptide_set) %>%
	group_by(qseqid,protein_acc, general_envo_name) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% lassopeptide_set) %>%
	group_by(general_envo_name) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

lasso_country <- gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% lassopeptide_set) %>%
	group_by(qseqid, protein_acc,geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

gmsc_samples_long_metadata_100 %>%
	filter(qseqid %in% lassopeptide_set) %>%
	filter(grepl("cow|cattle|rumen", general_envo_name, ignore.case = TRUE)) %>%
	group_by(qseqid, protein_acc, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

gmsc_samples_long_metadata_100 %>%
	filter(bg %in% "WP_002580678.1") %>%
	group_by(qseqid, protein_acc, general_envo_name, geographic_location) %>%
	summarise(n = n()) %>%
	arrange(desc(n))


########################################
# Plot world map of peptides at 100% ID
########################################
world <- map_data("world") 

unique(habitat_product$bgc_class)

map_plot <- ggplot() +
	geom_map(
		data = world, map = world,
		aes(x = long, y = lat, map_id = region),
		fill = "grey90", color = "black", size = 0.2
	) +
	geom_jitter(
		data = gmsc_samples_long_metadata_100 %>%
			filter(bgc_class %in% c(
				"Class II lanthipeptides",
				"Class I lanthipeptides",
				"Lassopeptides",
				"Core Peptides IIA",
				"Core Peptides IIC",
				"LAP",
				"Angicin",
				"Core Peptides IIB",
				"Thiopeptide",
				"Sactipeptide"
			)) %>%
			filter(grepl("rumen|gut", general_envo_name, ignore.case = TRUE)) %>%
			filter(!is.na(latitude) & !is.na(longitude)),
		aes(x = longitude, y = latitude, fill = bgc_class),
		size = 1,
		shape = 21,
		color = "black",
		alpha = 0.7,
		width = 1,
		height = 1
	) +
	scale_color_npg() +
	labs(x = "Longitude", y = "Latitude", color = "BGC class") +
	theme_classic() +
	labs(font = "Arial") +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.05, "cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  coord_fixed(
    1.3,
    xlim = c(-180, 180),
    ylim = c(-60, 85),  # adjust if you need more/less of the poles
    expand = FALSE
  ) +
  scale_fill_manual(
	values = c(
	  "Class II lanthipeptides" = rgb(154, 47, 45, maxColorValue = 255),
	  "Class I lanthipeptides" = rgb(239, 119, 0, maxColorValue = 255),
	  "Lassopeptides" = rgb(0, 60, 103, maxColorValue = 255),
	  "Core Peptides IIA" = rgb(54, 92, 141, maxColorValue = 255),
	  "Core Peptides IIC" = rgb(54, 92, 141, maxColorValue = 255),
	#   "LAP" = rgb(204, 102, 153, maxColorValue = 255),
	#   "Angicin" = rgb(54, 92, 141, maxColorValue = 255),
	  "Core Peptides IIB" = rgb(54, 92, 141,  maxColorValue = 255),
	  "Thiopeptide" = rgb(102, 51, 0, maxColorValue = 255),
	  "Sactipeptide" = rgb(0, 153, 76, maxColorValue = 255)
	),
	name = "BGC class"
  )


ggsave(
	"reports/figures/world_map_core_peptides_100ID.png",
	map_plot,
	width = 15,
	height = 10,
	units = "cm",
	dpi = 600
)


########################################
# Global summary stats on which peptides are in the most samples
########################################

gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class) %>%
	summarise(n = n()) %>%
	arrange(desc(n))

gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class, 
		# env
		microontology
		) %>%
	summarise(n = n()) %>%
	arrange(desc(n))  %>% View

gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class, 
		# env
		general_envo_name
		) %>%
	summarise(n = n()) %>%
	arrange(desc(n))  %>% View

gmsc_samples_long_metadata_100 %>%
	group_by(bgc_class, 
		# env
		general_envo_name
		) %>%
	summarise(n = n()) %>%
	arrange(desc(n))  %>% View


fwrite(gmsc_samples_long_metadata_100 %>%
	group_by(qseqid, bgc_class) %>%
	summarise(n = n()) %>%
	arrange(desc(n)),
	"reports/peptide_sample_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



########################################
# On average how many general_envo_name is each peptide in, by class?
########################################
gmsc_samples_long_metadata_100 %>%
	group_by(qseqid, bgc_class) %>%
	summarise(n_biomes = n_distinct(general_envo_name)) %>%
	ungroup() %>%
	group_by(bgc_class) %>%
	summarise(mean_biomes = mean(n_biomes), median_biomes = median(n_biomes), n_peptides = n()) %>%
	arrange(desc(mean_biomes))


gmsc_samples_long_metadata_100 %>%
	group_by(qseqid, bgc_class) %>%
	summarise(n_geo = n_distinct(geographic_location)) %>%
	ungroup() %>%
	group_by(bgc_class) %>%
	summarise(mean_geo = mean(n_geo), median_geo = median(n_geo), n_peptides = n()) %>%
	arrange(desc(mean_geo))
