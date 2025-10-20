#' @param: file_amide = table of Amide metabolites
#' @param: file_C8 = table of C8-pos metabolites
#' @param: file_C18 = table of C18-neg metabolites
#' @param: file_HILIC = table of HILIC metabolites



################################################################################
#.                      Things to look out for.                                #
#   #
################################################################################

build_mapping_file_function <- function(path_amide, path_C8, path_C18, path_HILIC)
  
{
  
  #--------------------------------------------------------------------------------#
  #---------------Read in files and select variables  -----------------------------#
  #--------------------------------------------------------------------------------#
  
  #########
  #Amide
  #########
  
  amide_info <- data.table::fread(path_amide) |>
    dplyr::select(Metabolite, DB_ID, Method, `Presence MESA X01`, `Presence MESA PILOT`, `Presence MESA2`) |> 
    dplyr::rename(Included_X01 = `Presence MESA X01`, 
                  Included_Pilot = `Presence MESA PILOT`, 
                  Included_MESA2 = `Presence MESA2`,
                  HMDB_ID = DB_ID) |>
    dplyr::mutate(Known_Compound = dplyr::case_when(!is.na(Metabolite) ~ 1, 
                                                    TRUE ~ 0),
                  Compound_ID_X01 = NA_character_) |>
    dplyr::mutate(Original_Metabolite_Name = Metabolite) |>
    dplyr::mutate(Metabolite = paste("Amide", make.names(Original_Metabolite_Name), sep="_")) |>
    dplyr::select(Metabolite, Compound_ID_X01, Original_Metabolite_Name, HMDB_ID, Known_Compound, Included_X01, Included_Pilot, Included_MESA2)
  
  amide_knowns <- as.data.frame(table(amide_info$Known_Compound))
  
  #########
  # C8
  #########
  
  C8_info <- data.table::fread(path_C8) |>
    dplyr::select(Metabolite, HMDB_ID, Method, Compound_ID_X01) |> 
    dplyr::mutate(Metabolite= dplyr::case_when(Compound_ID_X01=="TF04" ~ "LPC 16:0/0:0_v2", 
                                               TRUE ~ Metabolite)) |>
    dplyr::mutate(Included_X01 = 1, 
                  Included_Pilot = 1, 
                  Included_MESA2 = 1) 
  
  #Clean C8 metabolite names
  C8_info$Metabolite <- ifelse(is.na(C8_info$Metabolite), NA, ifelse(C8_info$Metabolite=="", NA, ifelse(C8_info$Metabolite=="NA", NA, C8_info$Metabolite)))
  
  #Continue formatting:
  C8_info <- C8_info |>
  dplyr::mutate(Known_Compound = dplyr::case_when(!is.na(Metabolite) ~ 1, 
                                                    TRUE ~ 0)) |>
    dplyr::mutate(Original_Metabolite_Name = Metabolite) |>
    dplyr::mutate(Metabolite = dplyr::case_when(!is.na(Metabolite) ~ paste("C8", make.names(Original_Metabolite_Name), sep="_"),
                                                TRUE ~ paste("C8",Compound_ID_X01, sep="_"))) |>
    dplyr::select(Metabolite, Original_Metabolite_Name, Compound_ID_X01, HMDB_ID, Known_Compound, Included_X01, Included_Pilot, Included_MESA2)
  
  #Knowns vs unkowns N
  C8_knowns <- as.data.frame(table(C8_info$Known_Compound))
  
  #########
  # C18
  #########
  
  C18_info <- data.table::fread(path_C18) |>
    dplyr::select(Metabolite, HMDB_ID, Method, Compound_ID_X01) |> 
    dplyr::mutate(Included_X01 = 1, 
                  Included_Pilot = 1, 
                  Included_MESA2 = 1) 
  
  #Clean C18 metabolite names
  C18_info$Metabolite <- ifelse(is.na(C18_info$Metabolite), NA, ifelse(C18_info$Metabolite=="", NA, ifelse(C18_info$Metabolite=="NA", NA, C18_info$Metabolite)))
  
  #Continue formatting:
  C18_info <- C18_info |>
    dplyr::mutate(Known_Compound = dplyr::case_when(!is.na(Metabolite) ~ 1, 
                                                    TRUE ~ 0)) |>
    dplyr::mutate(Original_Metabolite_Name = Metabolite) |>
    dplyr::mutate(Metabolite = dplyr::case_when(!is.na(Metabolite) ~ paste("C18", make.names(Original_Metabolite_Name), sep="_"),
                                                TRUE ~ paste("C18",Compound_ID_X01, sep="_"))) |>
    dplyr::select(Metabolite, Original_Metabolite_Name, Compound_ID_X01, HMDB_ID, Known_Compound, Included_X01, Included_Pilot, Included_MESA2)
  
  #Knowns vs unkowns N
  C18_knowns <- as.data.frame(table(C18_info$Known_Compound))
  
  #########
  # HILIC
  #########
  
  HILIC_info <- data.table::fread(path_HILIC) |>
    dplyr::select(Metabolite, HMDB_ID, Method, Compound_ID_X01) |> 
    dplyr::mutate(Included_X01 = 1, 
                  Included_Pilot = 1, 
                  Included_MESA2 = 1) 
  
  #Clean C18 metabolite names
  HILIC_info$Metabolite <- ifelse(is.na(HILIC_info$Metabolite), NA, ifelse(HILIC_info$Metabolite=="", NA, ifelse(HILIC_info$Metabolite=="NA", NA, HILIC_info$Metabolite)))
  
  #Continue formatting:
  HILIC_info <- HILIC_info |>
    dplyr::mutate(Known_Compound = dplyr::case_when(!is.na(Metabolite) ~ 1, 
                                                    TRUE ~ 0)) |>
    dplyr::mutate(Original_Metabolite_Name = Metabolite) |>
    dplyr::mutate(Metabolite = dplyr::case_when(!is.na(Metabolite) ~ paste("HILIC", make.names(Original_Metabolite_Name), sep="_"),
                                                TRUE ~ paste("HILIC",Compound_ID_X01, sep="_"))) |>
    dplyr::select(Metabolite, Original_Metabolite_Name, Compound_ID_X01, HMDB_ID, Known_Compound, Included_X01, Included_Pilot, Included_MESA2)
  
  #Knowns vs unkowns N
  HILIC_knowns <- as.data.frame(table(HILIC_info$Known_Compound))
  
  #########
  # Bind
  #########
  
  mapping <- rbind(amide_info, C8_info)
  mapping <- rbind(mapping, C18_info)
  mapping <- rbind(mapping, HILIC_info)
  
  all_knowns <- as.data.frame(table(mapping$Known_Compound))
  #----------------Build QC info------------------------#
  
  QC_info <-  list(
    knowns = list(amide_knowns = amide_knowns,
                  C8_knowns = C8_knowns,
                  C18_knowns = C18_knowns,
                  HILIC_knowns = HILIC_knowns,
                  all_knowns = all_knowns)
  )
  
  #----------------Outputs -----------------------------#
  
  list(
    
    QC_info_mapping_out = QC_info,
    mapping_final = mapping
    
  )
  
  
}
  
