#For testing:

#path_amide =  "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/25_0107_TOPMed_MESA_Amide-neg_rev031325.csv"
#path_C8 = "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C8-pos_checksums_rev031325.csv"
#path_C18 = "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C18-neg_checksums_rev031325.csv"
#path_HILIC= "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_HILIC-pos_checksums_rev031325.csv"

#path_amide_info = "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_AmideNeg_SampleInfo_20250329.txt"
#path_C8_info = "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C8Pos_SampleInfo_20250329.txt"
#path_C18_info = "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C18Neg_SampleInfo_20250329.txt"
#path_HILIC_info = "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_HILIC-Pos_SampleInfo_20250329.txt"
#sample_type_col = "Sample_type"
#QC_label = "QC-pooled_ref"
#####

QC_metabolites_function <- function(path_amide, path_C8, path_C18, path_HILIC, path_amide_info, path_C8_info, path_C18_info, path_HILIC_info, QC_label, Metabs_long)
  
{
  #Abundance tables
  amide_raw <- data.table::fread(path_amide)
  C8_raw <- data.table::fread(path_C8)
  C18_raw <- data.table::fread(path_C18)
  HILIC_raw <- data.table::fread(path_HILIC)
  
  #Info tables
  amide_info <- read.table(path_amide_info, header=TRUE, sep="\t")
  C8_info <- read.table(path_C8_info, header=TRUE, sep="\t")
  C18_info <- read.table(path_C18_info, header=TRUE, sep="\t")
  HILIC_info <- read.table(path_HILIC_info, header=TRUE, sep="\t")
  
  #----------------------------------------------------------#
  #---------------Amide QC info------------------------------#
  #----------------------------------------------------------#
  qc_cols_amide <- amide_info %>%
    filter(Sample_type == QC_label) %>%
    mutate(TOM_ID = trimws(TOM_ID)) %>%
    distinct(TOM_ID, .keep_all = FALSE) %>%
    filter(!is.na(TOM_ID) & TOM_ID != "") %>%
    pull(TOM_ID)
  
 message(length(qc_cols_amide), " pooled QC IDs identified: ") 
  
  amide_QC <- amide_raw |>
    dplyr::select(Metabolite, all_of(qc_cols_amide)) |>
    dplyr::mutate(key = paste("Amide", make.names(Metabolite), sep="_"),
                  Known = 1,
                  Assay = "Amide") |>
    dplyr::mutate(key <- make.names(key)) |>
    dplyr::mutate(Metabolite = key) |>
    dplyr::select(-key) |>
    dplyr::rowwise() %>%
    dplyr::mutate(
      mean_intensity = mean(c_across(all_of(qc_cols_amide)), na.rm = TRUE),
      sd_intensity   = sd(c_across(all_of(qc_cols_amide)), na.rm = TRUE),
      cv_percent     = (sd_intensity / mean_intensity) * 100,
      n_pools        = sum(!is.na(c_across(all_of(qc_cols_amide))))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Metabolite, Assay, Known, cv_percent, n_pools)
  
  #----------------------------------------------------------#
  #---------------C8 QC info------------------------------#
  #----------------------------------------------------------#
  
  qc_cols_C8 <- C8_info %>%
    filter(Sample_type == QC_label) %>%
    mutate(TOM_ID = trimws(TOM_ID)) %>%
    distinct(TOM_ID, .keep_all = FALSE) %>%
    filter(!is.na(TOM_ID) & TOM_ID != "") %>%
    pull(TOM_ID)
  
  message(length(qc_cols_C8), " pooled QC IDs identified: ") 
  
  C8_QC <- C8_raw |>
    dplyr::select(Metabolite, Compound_ID_X01, all_of(qc_cols_C8)) |>
    dplyr::mutate(
      Metabolite = dplyr::case_when(
        Compound_ID_X01 == "TF04" ~ "LPC 16:0/0:0_v2",
        TRUE ~ Metabolite
      ),
      Metabolite = dplyr::na_if(Metabolite, ""),
      Metabolite = dplyr::na_if(Metabolite, "NA"),
      Original_Metabolite_Name = Metabolite,
      Known = dplyr::if_else(is.na(Metabolite), 0, 1),
      Assay = "C8"
    ) |>
    dplyr::mutate(
      Metabolite = dplyr::if_else(
        !is.na(Metabolite),
        paste0("C8_", make.names(Original_Metabolite_Name)),
        paste0("C8_", Compound_ID_X01)
      )
    ) |>
    dplyr::mutate(
      mean_intensity = rowMeans(dplyr::pick(all_of(qc_cols_C8)), na.rm = TRUE),
      sd_intensity   = apply(
        as.data.frame(dplyr::pick(all_of(qc_cols_C8))),
        1,
        sd,
        na.rm = TRUE
      ),
      cv_percent     = (sd_intensity / mean_intensity) * 100,
      n_pools        = rowSums(!is.na(dplyr::pick(all_of(qc_cols_C8))))
    ) |>
    dplyr::select(Metabolite, Assay, Known, cv_percent, n_pools)
  
  
  

  #----------------------------------------------------------#
  #---------------C18 QC info------------------------------#
  #----------------------------------------------------------#
  
  qc_cols_C18 <- C18_info %>%
    filter(Sample_type == QC_label) %>%
    mutate(TOM_ID = trimws(TOM_ID)) %>%
    distinct(TOM_ID, .keep_all = FALSE) %>%
    filter(!is.na(TOM_ID) & TOM_ID != "") %>%
    pull(TOM_ID)
  
  message(length(qc_cols_C18), " pooled QC IDs identified: ") 
  
  C18_QC <- C18_raw |>
    dplyr::select(Metabolite, Compound_ID_X01, all_of(qc_cols_C18)) |>
    dplyr::mutate(Metabolite = dplyr::case_when(is.na(Metabolite) ~ NA_character_,
                                                Metabolite=="" ~ NA_character_,
                                                Metabolite=="NA" ~ NA_character_,
                                                TRUE ~ Metabolite)) |>
    dplyr::mutate(Original_Metabolite_Name = Metabolite,
                  Known = dplyr::case_when(is.na(Metabolite) ~ 0,
                                           Metabolite=="" ~ 0,
                                           Metabolite=="NA" ~ 0,
                                           TRUE ~ 1),
                  ,
                  Assay = "C18") |>
    dplyr::mutate(Metabolite = dplyr::case_when(!is.na(Metabolite) ~ paste("C18", make.names(Original_Metabolite_Name), sep="_"),
                                                TRUE ~ paste("C18",Compound_ID_X01, sep="_")))|>
    dplyr::rowwise() %>%
    dplyr::mutate(mean_intensity = mean(c_across(all_of(qc_cols_C18)), na.rm = TRUE),
                  sd_intensity   = sd(c_across(all_of(qc_cols_C18)), na.rm = TRUE),
                  cv_percent     = (sd_intensity / mean_intensity) * 100,
                  n_pools        = sum(!is.na(c_across(all_of(qc_cols_C18))))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Metabolite, Assay, Known,  cv_percent, n_pools)
  
  
  #----------------------------------------------------------#
  #---------------HILIC QC info------------------------------#
  #----------------------------------------------------------#
  qc_cols_HILIC <- HILIC_info %>%
    filter(Sample_type == QC_label) %>%
    mutate(TOM_ID = trimws(TOM_ID)) %>%
    distinct(TOM_ID, .keep_all = FALSE) %>%
    filter(!is.na(TOM_ID) & TOM_ID != "") %>%
    pull(TOM_ID)
  
  message(length(qc_cols_HILIC), " pooled QC IDs identified: ") 
  
  HILIC_QC <- HILIC_raw |>
    dplyr::select(Metabolite, Compound_ID_X01, all_of(qc_cols_HILIC)) |>
    dplyr::mutate(Metabolite = dplyr::case_when(is.na(Metabolite) ~ NA_character_,
                                                Metabolite=="" ~ NA_character_,
                                                Metabolite=="NA" ~ NA_character_,
                                                TRUE ~ Metabolite)) |>
    dplyr::mutate(Original_Metabolite_Name = Metabolite) |>
    dplyr::mutate(Metabolite = dplyr::case_when(!is.na(Metabolite) ~ paste("HILIC", make.names(Original_Metabolite_Name), sep="_"),
                                                TRUE ~ paste("HILIC",Compound_ID_X01, sep="_")),
                  Known = dplyr::case_when(is.na(Metabolite) ~ 0,
                           Metabolite=="" ~ 0,
                           Metabolite=="NA" ~ 0,
                           TRUE ~ 1),
                  Assay = "HILIC") |>
    
    dplyr::rowwise() %>%
    dplyr::mutate(mean_intensity = mean(c_across(all_of(qc_cols_HILIC)), na.rm = TRUE),
                  sd_intensity   = sd(c_across(all_of(qc_cols_HILIC)), na.rm = TRUE),
                  cv_percent     = (sd_intensity / mean_intensity) * 100,
                  n_pools        = sum(!is.na(c_across(all_of(qc_cols_HILIC))))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Metabolite, Assay, Known,  cv_percent, n_pools)
  

  
  #----------------------------------------------------------#
  #---------------final_file---------------------------------#
  #----------------------------------------------------------#
  
  CV_file <- rbind(amide_QC, C8_QC)
  CV_file <- rbind(CV_file, C18_QC)
  CV_file <- rbind(CV_file, HILIC_QC)
  
  #-------------------------------------------------------------#
  #---------------Missing info---------------------------------#
  #------------------------------------------------------------#
  
  metab_cols <- Metabs_long |>
    dplyr::select(-idno, -sidno, -TOM_ID, -subject_id, -exam) |>
    names()
  

    
    ##Long version:
    
    Miss_long <- Metabs_long |>
      dplyr::group_by(exam) |>
      dplyr::summarise(across(all_of(metab_cols), ~ 100 * mean(is.na(.)), .names = "{.col}")) |>
      dplyr::bind_rows(
        Metabs_long |>
          dplyr::summarise(across(all_of(metab_cols), ~ 100 * mean(is.na(.)), .names = "{.col}")) |>
          dplyr::mutate(exam = 0)
      ) |>
      tidyr::pivot_longer(
        cols = all_of(metab_cols),
        names_to = "Metabolite",
        values_to = "percent_missing"
      ) |>
      tidyr::pivot_wider(names_from = exam, names_prefix="missing_exam_", values_from = percent_missing) |>
      dplyr::rename(missing_all_exams = missing_exam_0)

    #--------------------Build----------------#
    
    final_file <- dplyr::full_join(CV_file,Miss_long, dplyr::join_by(Metabolite))
    
    #--------------------Outputs----------------#


  list (
QC_info = final_file
        )

  
}