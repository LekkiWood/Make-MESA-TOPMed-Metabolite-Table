#' @param: file_amide = table of Amide metabolites
#' @param: file_C8 = table of C8-pos metabolites
#' @param: file_C18 = table of C18-neg metabolites
#' @param: file_HILIC = table of HILIC metabolites
#' @param: file_bridge = Bridging file
#' @param: path_amide_info = file with sampleorder info, TOPMed ID, exam, injection order, sample type, raw file name, project and date_injected for amide assay.
#' @param: path_C8_info = file with sampleorder info, TOPMed ID, exam, injection order, sample type, raw file name, project and date_injected for C8-pos assay.
#' @param: path_C18_info = file with sampleorder info, TOPMed ID, exam, injection order, sample type, raw file name, project and date_injected for C18-neg assay.
#' @param: path_HILIC_info = file with sampleorder info, TOPMed ID, exam, injection order, sample type, raw file name, project and date_injected for HILIC assay.


################################################################################
#.                      Things to look out for.                                #
# Check the number of delimiters in the quarto file for reading in filenames.  #
################################################################################


build_metabs_table_function <- function(path_amide, path_C8, path_C18, path_HILIC, path_amide_info, path_C8_info, path_C18_info, path_HILIC_info, path_bridge)
  
{
  
  #----------------------------------------------------------#
  #---------------Read in files -----------------------------#
  #----------------------------------------------------------#
  
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
  
  #Bridging file
  bridge <- data.table::fread(path_bridge) |>
    dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
    dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`)
  
  #Filenames saved below for calling to Quarto file
  
  #----------------------------------------------------------#
  #---------------Clean Amide -------------------------------#
  #----------------------------------------------------------#
  
  amide_info_TOMIDs <- amide_info |>
    dplyr::filter(substr(TOM_ID,1,3) == "TOM")
  
  sample_info_dims_amide <- dim(amide_info_TOMIDs)
  
  amide <- amide_raw |>
    dplyr::select(Metabolite, all_of(amide_info_TOMIDs$TOM_ID))
  
  column_names_amide = names(amide)
  duplicate_indices_amide <- duplicated(column_names_amide)
  duplicate_column_names_amide <- column_names_amide[duplicate_indices_amide]
  
  #No duplicates
  
  phase1_dims_amide <- dim(amide)
  
  amide$TOM148122 <- ifelse(amide$TOM148122=="1+2559:2584.03935869994933", NA, amide$TOM148122)
  amide$TOM148122 <- as.numeric(amide$TOM148122)
  
  
  amide$key <- paste("Amide", make.names(amide$Metabolite), sep="_")
  amide$key <- make.names(amide$key)
  amide$Metabolite<- amide$key
  amide$key <- NULL
  
  
  #check distributions of vars 
  vartypes_amide <- as.data.frame(t(as.data.frame(lapply(amide, class))))
  names(vartypes_amide) <- "Freq"
  types_amide <- as.data.frame(table(vartypes_amide$Freq))
  #6 are logical due to missing data
  
  
  amide <- amide |>
    tidyr::pivot_longer(cols=c(-Metabolite),names_to="TOM_ID") |>
    tidyr::pivot_wider(names_from=c(Metabolite)) 
  
  
  amide_merged <- amide_info_TOMIDs |>
    dplyr::select(TOM_ID, exam, subject_id) |> 
    dplyr::mutate(sidno = subject_id) |>
    dplyr::full_join(amide, dplyr::join_by(TOM_ID))
  
  
  #Identify duplicates
  keys <- c("sidno", "exam")
  # Which side(s) has non-unique key combos?
  amide_dups <- amide_merged |>
    dplyr::count(across(all_of(keys)), name = "n") |> 
    dplyr::filter(n > 1)
  
  amide_final <- amide_merged |>
    dplyr::filter(dplyr::n() == 1, .by = c(sidno, exam))
  
  final_table_dims_amide <- dim(amide_final)
  
  #----------------------------------------------------------#
  #------------------Clean C8 -------------------------------#
  #----------------------------------------------------------#
  
  C8_info_TOMIDs <- C8_info |>
    dplyr::filter(substr(TOM_ID,1,3) == "TOM") 
  
  sample_info_dims_C8 <- dim(C8_info_TOMIDs)
  
  C8 <- C8_raw |>
    dplyr::select(Metabolite, Compound_ID_X01, all_of(C8_info_TOMIDs$TOM_ID)) |>
    dplyr::mutate(Metabolite= dplyr::case_when(Compound_ID_X01=="TF04" ~ "LPC 16:0/0:0_v2", 
                                               TRUE ~ Metabolite))
  
  #check_integer <- as.data.frame(C8[,names(C8)=="Metabolite" | sapply(C8,is.integer64)])
  #C8_integer <- as.data.frame(C8[, sapply(C8,is.integer64)])
  
  #for (i in names(C8_integer)){
  
  
  #  C8[[i]] <- as.double(C8[[i]])
  #}
  
  column_names_C8 = names(C8)
  duplicate_indices_C8 <- duplicated(column_names_C8)
  duplicate_column_names_C8 <- column_names_C8[duplicate_indices_C8]
  
  #No duplicates
  
  phase1_dims_C8 <- dim(C8)
  
  C8$Metabolite <- ifelse(is.na(C8$Metabolite), NA, ifelse(C8$Metabolite=="", NA, ifelse(C8$Metabolite=="NA", NA, C8$Metabolite)))
  C8$key <- ifelse(is.na(C8$Metabolite), paste("C8", make.names(C8$Compound_ID_X01), sep="_"), paste("C8", make.names(C8$Metabolite), sep="_"))
  
  #C8$Metabolite <- ifelse(C8$Metabolite=="", NA,  C8$Metabolite)
  #C8$Metabolite <- make.names(C8$Metabolite, unique = TRUE)
  #C8$key <- ifelse(is.na(C8$Metabolite), paste("C8", make.names(C8$Compound_ID_X01), sep="_"), paste("C8", make.names(C8$Metabolite), sep="_"))
  #C8$key <- make.names(C8$key)
  C8$Metabolite<- C8$key
  C8$Metabolite<- make.names(C8$Metabolite)
  C8$key <- NULL
  
  #check distributions of vars 
  vartypes_C8 <- as.data.frame(t(as.data.frame(lapply(C8, class))))
  names(vartypes_C8) <- "Freq"
  types_C8 <- as.data.frame(table(vartypes_C8$Freq))
  
  C8 <- subset(C8, select= -Compound_ID_X01)
  
  
  C8 <-  data.table::transpose(C8, keep.names="TOM_ID", make.names="Metabolite")
  
  names(C8) <- make.names(names(C8), unique=TRUE)
  
  C8_merged <- C8_info_TOMIDs |>
    dplyr::select(TOM_ID, exam, subject_id) |> 
    dplyr::mutate(sidno = subject_id) |>
    dplyr::full_join(C8, dplyr::join_by(TOM_ID))
  
  
  
  
  #Identify duplicates
  keys <- c("sidno", "exam")
  # Which side(s) has non-unique key combos?
  C8_dups <- C8_merged |>
    dplyr::count(across(all_of(keys)), name = "n") |> 
    dplyr::filter(n > 1)
  
  C8_final <- C8_merged |>
    dplyr::filter(dplyr::n() == 1, .by = c(sidno, exam))
  
  final_table_dims_C8 <- dim(C8_final)  
  
  #----------------------------------------------------------#
  #------------------Clean C18 -------------------------------#
  #----------------------------------------------------------#
  
  C18_info_TOMIDs <- C18_info |>
    dplyr::filter(substr(TOM_ID,1,3) == "TOM")
  
  sample_info_dims_C18 <- dim(C18_info_TOMIDs)
  
  C18 <- C18_raw |>
    dplyr::select(Metabolite, Compound_ID_X01, all_of(C18_info_TOMIDs$TOM_ID))
  
  
  column_names_C18 = names(C18)
  duplicate_indices_C18 <- duplicated(column_names_C18)
  duplicate_column_names_C18 <- column_names_C18[duplicate_indices_C18]
  
  #No duplicates
  
  phase1_dims_C18 <- dim(C18)
  
  C18$Metabolite <- ifelse(is.na(C18$Metabolite), NA, ifelse(C18$Metabolite=="", NA, ifelse(C18$Metabolite=="NA", NA, C18$Metabolite)))
  C18$key <- ifelse(is.na(C18$Metabolite), paste("C18", make.names(C18$Compound_ID_X01), sep="_"), paste("C18", make.names(C18$Metabolite), sep="_"))
  
  C18$Metabolite<- C18$key
  C18$Metabolite<- make.names(C18$Metabolite)
  C18$key <- NULL
  
  #check distributions of vars 
  vartypes_C18 <- as.data.frame(t(as.data.frame(lapply(C18, class))))
  names(vartypes_C18) <- "Freq"
  types_C18 <- as.data.frame(table(vartypes_C18$Freq))
  
  C18 <- subset(C18, select= -Compound_ID_X01)
  
  
  C18 <-  data.table::transpose(C18, keep.names="TOM_ID", make.names="Metabolite")
  
  names(C18) <- make.names(names(C18), unique=TRUE)
  
  C18_merged <- C18_info_TOMIDs |>
    dplyr::select(TOM_ID, exam, subject_id) |> 
    dplyr::mutate(sidno = subject_id) |>
    dplyr::full_join(C18, dplyr::join_by(TOM_ID))
  
  
  
  #Identify duplicates
  keys <- c("sidno", "exam")
  # Which side(s) has non-unique key combos?
  C18_dups <- C18_merged |>
    dplyr::count(across(all_of(keys)), name = "n") |> 
    dplyr::filter(n > 1)
  
  C18_final <- C18_merged |>
    dplyr::filter(dplyr::n() == 1, .by = c(sidno, exam))
  
  final_table_dims_C18 <- dim(C18_final)  
  
  #-------------------------------------------------------------#
  #------------------Clean HILIC -------------------------------#
  #-------------------------------------------------------------#
  
  HILIC_info_TOMIDs <- HILIC_info |>
    dplyr::filter(substr(TOM_ID,1,3) == "TOM")
  
  sample_info_dims_HILIC <- dim(HILIC_info_TOMIDs)
  
  HILIC <- HILIC_raw |>
    dplyr::select(Metabolite, Compound_ID_X01, all_of(HILIC_info_TOMIDs$TOM_ID))
  
  #check_integer <- as.data.frame(HILIC[,names(HILIC)=="Metabolite" | sapply(HILIC,is.integer64)])
  #HILIC_integer <- as.data.frame(HILIC[, sapply(HILIC,is.integer64)])
  
  #for (i in names(HILIC_integer)){
  
  
  #  HILIC[[i]] <- as.double(HILIC[[i]])
  #}
  
  column_names_HILIC = names(HILIC)
  duplicate_indices_HILIC <- duplicated(column_names_HILIC)
  duplicate_column_names_HILIC <- column_names_HILIC[duplicate_indices_HILIC]
  
  #No duplicates
  
  phase1_dims_HILIC <- dim(HILIC)
  
  HILIC$Metabolite <- ifelse(is.na(HILIC$Metabolite), NA, ifelse(HILIC$Metabolite=="", NA, ifelse(HILIC$Metabolite=="NA", NA, HILIC$Metabolite)))
  HILIC$key <- ifelse(is.na(HILIC$Metabolite), paste("HILIC", make.names(HILIC$Compound_ID_X01), sep="_"), paste("HILIC", make.names(HILIC$Metabolite), sep="_"))
  
  #HILIC$Metabolite <- ifelse(HILIC$Metabolite=="", NA,  HILIC$Metabolite)
  #HILIC$Metabolite <- make.names(HILIC$Metabolite, unique = TRUE)
  #HILIC$key <- ifelse(is.na(HILIC$Metabolite), paste("HILIC", make.names(HILIC$Compound_ID_X01), sep="_"), paste("HILIC", make.names(HILIC$Metabolite), sep="_"))
  #HILIC$key <- make.names(HILIC$key)
  HILIC$Metabolite<- HILIC$key
  HILIC$Metabolite<- make.names(HILIC$Metabolite)
  HILIC$key <- NULL
  
  #check distributions of vars 
  vartypes_HILIC <- as.data.frame(t(as.data.frame(lapply(HILIC, class))))
  names(vartypes_HILIC) <- "Freq"
  types_HILIC <- as.data.frame(table(vartypes_HILIC$Freq))
  
  HILIC <- subset(HILIC, select= -Compound_ID_X01)
  
  
  HILIC <-  data.table::transpose(HILIC, keep.names="TOM_ID", make.names="Metabolite")
  
  names(HILIC) <- make.names(names(HILIC), unique=TRUE)
  
  HILIC_merged <- HILIC_info_TOMIDs |>
    dplyr::select(TOM_ID, exam, subject_id) |> 
    dplyr::mutate(sidno = subject_id) |>
    dplyr::full_join(HILIC, dplyr::join_by(TOM_ID))
  
  
  
  
  #Identify duplicates
  keys <- c("sidno", "exam")
  # Which side(s) has non-unique key combos?
  HILIC_dups <- HILIC_merged |>
    dplyr::count(across(all_of(keys)), name = "n") |> 
    dplyr::filter(n > 1)
  
  HILIC_final <- HILIC_merged |>
    dplyr::filter(dplyr::n() == 1, .by = c(sidno, exam))
  
  final_table_dims_HILIC <- dim(HILIC_final)  
  
  #----------------Merge-----------------------------#
  
  
  metabs_final <- merge(amide_final, C8_final, by=c("sidno", "exam", "subject_id", "TOM_ID"), all=TRUE)
  metabs_final <- merge(metabs_final, C18_final, by=c("sidno", "exam", "subject_id", "TOM_ID"), all=TRUE)
  metabs_final <- merge(metabs_final, HILIC_final, by=c("sidno", "exam", "subject_id", "TOM_ID"), all=TRUE)
  
  metabs_final_dim <- dim(metabs_final)
  metabs_final_uniqueid <- length(unique(metabs_final$sidno))
  metabs_final_uniqueid_E1 <- length(unique(metabs_final[which(metabs_final$exam==1),]$sidno))
  metabs_final_uniqueid_E5 <- length(unique(metabs_final[which(metabs_final$exam==5),]$sidno))
  metabs_final_uniqueid_E6 <- length(unique(metabs_final[which(metabs_final$exam==6),]$sidno))
  
  miss_sidno1 <- metabs_final$sidno
  miss_sidno <- miss_sidno1[!miss_sidno1 %in% bridge$sidno]
  
  metabs_final <- dplyr::left_join(metabs_final, bridge, dplyr::join_by(sidno))
  final_metabs_merged_dim <- dim(metabs_final)
  
  
  #----------------Build QC info------------------------#
  
  QC_info <- list(
    filenames = list(raw_amide_file_name = path_amide,
                     raw_C8_file_name = path_C8,
                     raw_C18_file_name = path_C18,
                     raw_HILIC_file_name = path_HILIC,
                     
                     raw_amide_info_file_name = path_amide_info,
                     raw_C8_info_file_name = path_C8_info,
                     raw_C18_info_file_name = path_C18_info,
                     raw_HILIC_info_file_name = path_HILIC_info,
                     
                     raw_bridgingfile_file_name = path_bridge
                     
    ),
    
    sample_info_dims = list(
      sample_info_dims_amide = sample_info_dims_amide,
      sample_info_dims_C8 = sample_info_dims_C8,
      sample_info_dims_C18 = sample_info_dims_C18,
      sample_info_dims_HILIC = sample_info_dims_HILIC
    ),
    
    phase1_dims = list(
      amide_phase1_dims = phase1_dims_amide,
      C8_phase1_dims = phase1_dims_C8,
      C18_phase1_dims = phase1_dims_C18,
      HILIC_phase1_dims = phase1_dims_HILIC
    ),
    
    duplicate_TOMIDs = list(
      duplicate_column_names_amide = duplicate_column_names_amide,
      duplicate_column_names_C8 = duplicate_column_names_C8,
      duplicate_column_names_C18 = duplicate_column_names_C18,
      duplicate_column_names_HILIC = duplicate_column_names_HILIC
    ),
    
    duplicates_in_raw_tables = list(
      amide_dupes = amide_dups,
      C8_dupes = C8_dups,
      C18_dupes = C18_dups,
      HILIC_dupes = HILIC_dups
    ),
    
    final_table_dims = list(
      final_table_dims_amide= final_table_dims_amide,
      final_table_dims_C8= final_table_dims_C8,
      final_table_dims_C18= final_table_dims_C18,
      final_table_dims_HILIC= final_table_dims_HILIC
    ),
    final_df_info = list(
      metabs_final_dim = metabs_final_dim,
      metabs_final_uniqueid = metabs_final_uniqueid,
      metabs_final_uniqueid_E1= metabs_final_uniqueid_E1 ,
      metabs_final_uniqueid_E5  = metabs_final_uniqueid_E5,
      metabs_final_uniqueid_E6  = metabs_final_uniqueid_E6,
      miss_sidno = miss_sidno,
      final_metabs_merged_dim = final_metabs_merged_dim,
      final_metabs_sidno = miss_sidno1
    )
  )
  
  
  
  #----------------Outputs -----------------------------#
  
  list(
    
    QC_info_out = QC_info,
    metabs_final = metabs_final
    
  )
  
  
}