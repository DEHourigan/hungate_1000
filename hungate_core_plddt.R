# bacteriocins comparison
# Load libs
  .libPaths()
  .libPaths( c( "/data/san/data0/users/david/rstudio/packages" , .libPaths() ) )
  newlib <- "/data/san/data0/users/david/rstudio/packages"
  .libPaths()

  #install.packages("bio3d")
  library(bio3d)
  library(ggplot2)
  library(gplots)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  library(data.table)  
  library(dplyr)
  library("rjson")
  library(tidyr)


setwd("/data/san/data1/users/david/hungate_proteins/alphafold")
bacteriocin_dir = "/data/san/data1/users/david/hungate_proteins/alphafold"
  # Get all directories and subdirectories
bacteriocin_dirs <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE,
    pattern = "debug.json")

  # Initialize a list to store the data
  data_list_b <- list()

  # Loop over each file path and read the JSON file
  for(json_file in bacteriocin_dirs) {
    if(file.exists(json_file)){
      json_data <- fromJSON(file=json_file) %>% as.data.frame()

      # Create a unique id for each json file, based on its path
      file_id <- gsub(paste0("^", gsub("/", "\\/", bacteriocin_dir), "\\/"), "", dirname(json_file))
      file_id <- gsub("/", "_", file_id)  # Replace any remaining slashes with underscores

      # Add the id column to the dataframe
      json_data$id <- file_id

      # Store the data in the list, using the file_id as the key
      data_list_b[[file_id]] <- json_data
    } else {
      warning(paste("File does not exist:", json_file))
    }
  }


  # Combine all data frames into a single one
  final_data_bacteriocin <- bind_rows(data_list_b)
  bacteriocin_df_long <- final_data_bacteriocin %>%  pivot_longer(cols = starts_with("plddts"),
                 names_to = "prediction",
                 values_to = "Value")

  bacteriocin_df_long_clean = bacteriocin_df_long %>% group_by(id) %>%
                          mutate(mean_plddt = mean(Value)) %>% 
                          select(-order) %>% 
                          distinct() 


  # Assuming bacteriocin_df_long_clean is your data frame and it has a column named 'id'
  #bacteriocin_df_long_clean$id <- sub("^([A-Za-z0-9]+_[A-Za-z0-9]+)_.*", "\\1", bacteriocin_df_long_clean$id)

  # Adjust the pattern in sub to handle cases with a suffix like '-2'
  bacteriocin_df_long_clean$id_2 <- sub("^([A-Za-z0-9]+_[A-Za-z0-9]+(-[0-9]+)?).*", "\\1", bacteriocin_df_long_clean$id)
  bacteriocin_df_long_clean$id_3 <- gsub("circular_bacteriocin_|pediocin_like_", "", bacteriocin_df_long_clean$id)




  color_palette <- scales::col_numeric(c("tomato", "#FFFF00", "dodgerblue","blue"), domain = c(100, 70, 50))

  circular_core = c(
    "CDPAON_01140", # bingo
    "CECIIA_04450", # bingo
    #"DPGIHA_03245",   not a cluster
    "EAGELN_15015",  # clopol DSM
    "EBIEFA_02410", # bingo 
    "ENIMMA_01450", # bingo 
   # "GACLNI_20115",
    "GAPOLK_04680", # bingo
    "MFGIFH_02585", # bingo
    "NNIJIM_07585",  # bingo
    "NPEPEM_07125", # bingo 
    "POHBGK_08970", # bingo
    "GILPNI_10090",
    "LHAIJF_07965",
    "NLCFLM_27635"
  )
  bacteriocin_plddt_plot = bacteriocin_df_long_clean %>%
    #filter(id_3  %in% circular_core) %>%
    filter(!grepl("LHAIJF", id_3)) %>%
    filter(!grepl("class_IId", id_3)) %>%
    ggplot(aes(x = reorder(id_3, -mean_plddt), y = Value)) +
    geom_boxplot(aes(fill = id_2), outlier.shape = NA) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(
      values = c("#A73030FF", "#EFC000FF", "red","white"),
      guide = "none"
    ) +
    labs(y = "mean pLDDT", x = "Structural confidence") 


  ggsave(bacteriocin_plddt_plot, file = "plddt.png", width = 10, height = 10, units = "cm")
