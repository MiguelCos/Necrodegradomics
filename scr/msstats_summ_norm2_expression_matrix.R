# Function to get MSstats summarized and normalized output into wide format (expression matrix)  
# v 0.1 Miguel Cosenza  

library(dplyr)
library(tidyr)

msstats_out2_wide <- function(msstats_out) {
          
          tab_proccessed_data <- msstats_out$RunlevelData
          
          boxp1 <- ggplot(data = tab_proccessed_data,
                          aes(x = SUBJECT_ORIGINAL, y = LogIntensities)) +
                    geom_boxplot() + 
                    ggtitle("Log Protein Intensities by Sample from 'RunLevelData'")
          
          # select interesting columns to prep wide format table ----
          
          tab_selection_msts_data <- dplyr::select(tab_proccessed_data,
                                                   Protein, 
                                                   LogIntensities, 
                                                   Run = originalRUN, 
                                                   Group = GROUP_ORIGINAL,
                                                   BioReplicate = SUBJECT_ORIGINAL) %>%
                    mutate(GROUP_REPLICATE_RUN = paste0(Group,"_",BioReplicate)) %>%
                    dplyr::select(Protein, 
                                  LogIntensities, 
                                  GROUP_REPLICATE_RUN)
          
          # wide formating ----
          
          tab_wide_msts_data <- tidyr::pivot_wider(data = tab_selection_msts_data,
                                                   names_from = GROUP_REPLICATE_RUN, 
                                                   values_from = LogIntensities) %>% 
                    dplyr::rename(ID = Protein) %>%
                    dplyr::select(ID,
                                  matches("06h"),
                                  matches("12h"),
                                  matches("18h"),
                                  matches("24h"),
                                  matches("48h"),
                                  matches("72h"),
                                  matches("96h"))
          
}