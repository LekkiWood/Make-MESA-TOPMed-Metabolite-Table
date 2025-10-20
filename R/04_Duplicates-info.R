flag_duplicates_function <- function(QC_file, mapping_file)
  
{
    

###Old code to identify possible errors where strong does not match HMDB_ID
#errors <- Duplicates_file %>%
#  group_by(HMDB_ID) %>%
#  mutate(
    # 1️⃣ Flag same HMDB_ID but different names
#    Duplicate_flag = if_else(
#      n() > 1 & n_distinct(Cleaned_Metabolite_Name) > 1,
#      TRUE,
#      FALSE
#    ))


Duplicates_file <- mapping_file |>
  dplyr::select(Metabolite, Original_Metabolite_Name, HMDB_ID) |>
  dplyr::full_join(QC_file |> 
                     dplyr::select(Metabolite, cv_percent, Assay), dplyr::join_by(Metabolite)) |>
  dplyr::mutate(HMDB_ID = dplyr::case_when(Assay=="Amide" & Original_Metabolite_Name=="Glutamic acid" ~ "HMDB0000148",
                                           TRUE ~ HMDB_ID)) |>
  dplyr::mutate(Metabolite = dplyr::case_when(Assay=="Amide" & HMDB_ID=="HMDB0000122" ~ "Amide_Hexose",
                                             TRUE ~ Metabolite)) |>
  dplyr::mutate(cv_percent_temp = dplyr::case_when(is.na(cv_percent) ~ 100000,
                                                   TRUE ~ cv_percent)) |>
  dplyr::group_by(HMDB_ID) |>
  dplyr::mutate(Retain = dplyr::case_when(is.na(HMDB_ID) | HMDB_ID=="" ~ 0,
                                          n() == 1 ~ 0,
                                          cv_percent_temp== min(cv_percent_temp, na.rm = TRUE) ~ 1,
                                          HMDB_ID=="internal standard" ~ 2,
                                          # unique HMDB_ID      # lowest CV among duplicates
                                          TRUE ~ 3 ) )  |>
  dplyr::ungroup() |>
  dplyr::mutate(Retain = factor(Retain, labels=c("Unique or missing HMDB ID", "Duplicated HMDB ID with lowest CV", "Internal standard", "Duplicated HMDB_ID and not lowest CV"))) |>
  dplyr::select(-cv_percent_temp)

#table(Duplicates_file$Retain)

#------------------------Output----------------------#

#Duplicates_file

}

