library('targets')
library('tarchetypes')
library(crew)
tar_option_set(controller = crew_controller_local(workers = 10)) 
Sys.setenv(VROOM_CONNECTION_SIZE = as.character(10 * 1024 * 1024)) #For the large metab datasets

####This runs a workflow for cleaning the MESA TOPMed multi-omics project metabolite data, formatting it into the following format: wide for metabolites & long for exams, and providing basic quality control (QC) metrics



tar_option_set(packages = c("dplyr", "tidyr", "tibble", "readr", "data.table", "bit64", "foreign", "quarto", "rlang", "purrr", "rcompanion"))

tar_source("/media/Analyses/Make-MESA-TOPMed-Metabolite-Table/R")

list(
#---------------------------------------------------------------------------------------#
#--------------------------------1. Build metabolite table------------------------------#
#---------------------------------------------------------------------------------------#

#Abundance tables
tar_target(path_amide, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/25_0107_TOPMed_MESA_Amide-neg_rev031325.csv", format = "file"),
tar_target(path_C8, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C8-pos_checksums_rev031325.csv", format = "file"),
tar_target(path_C18, "/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_C18-neg_checksums_rev031325.csv", format = "file"),
tar_target(path_HILIC,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/24_1210_TOPMed_MESA_HILIC-pos_checksums_rev031325.csv", format = "file"),

#Info tables
tar_target(path_amide_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_AmideNeg_SampleInfo_20250329.txt", format = "file"),
tar_target(path_C8_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C8Pos_SampleInfo_20250329.txt", format = "file"),
tar_target(path_C18_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_C18Neg_SampleInfo_20250329.txt", format = "file"),
tar_target(path_HILIC_info,"/media/RawData/MESA/MESA-Multiomics/MESA-Multiomics_Metabolomics/MesaMetabolomics_PilotX01_HILIC-Pos_SampleInfo_20250329.txt", format = "file"),

#Bridging file
tar_target(path_bridge,"/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv", format = "file"),

#Run function
tar_target(build_metabs_out,   
           build_metabs_table_function(
             path_amide = path_amide,
             path_C8 = path_C8,
             path_C18 = path_C18,
             path_HILIC = path_HILIC,
             path_amide_info = path_amide_info,
             path_C8_info = path_C8_info,
             path_C18_info = path_C18_info,
             path_HILIC_info = path_HILIC_info,
             path_bridge = path_bridge)
),

#Save outputs
tar_target(build_metabs_QC, build_metabs_out$QC_info),
tar_target(Metabs_long, build_metabs_out$metabs_final),

#Save metabolite table (long form) as csv file
tar_target(long_metabtable_filename, paste0("MESA_TOPMed_Metabolite_Longform_", Sys.Date(), ".csv")),
tar_target(save_long_metab_table_csv,
  {
    out_dir  <- "outputs"
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(out_dir, long_metabtable_filename)
    readr::write_csv(Metabs_long, out_path)
    out_path
  },
  format = "file"
),

#---------------------------------------------------------------------------------------#
#--------------------------------2. Build mapping file.   ------------------------------#
#---------------------------------------------------------------------------------------#


tar_target(build_mapping_file_out, build_mapping_file_function(
  path_amide = path_amide,
  path_C8 = path_C8,
  path_C18 = path_C18,
  path_HILIC = path_HILIC)),

#Save outputs
tar_target(build_mapping_QC, build_mapping_file_out$QC_info_mapping_out),
tar_target(Mapping_file, build_mapping_file_out$mapping_final),

#Save mapping file as csv file
tar_target(mapping_filename, paste0("MESA_TOPMed_Metabolite_Mappingfile_", Sys.Date(), ".csv")),
tar_target(mapping_file_csv,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, mapping_filename)
             readr::write_csv(Mapping_file, out_path)
             out_path
           },
           format = "file"
),

#

#----------------------------------------------------------------------------#
#--------------------------------3. CV and Missingness-----------------------#
#----------------------------------------------------------------------------#


tar_target(make_metabolite_QC, QC_metabolites_function(path_amide = path_amide, 
                                                      path_C8 = path_C8, 
                                                      path_C18 = path_C18, 
                                                      path_HILIC = path_HILIC, 
                                                      path_amide_info = path_amide_info, 
                                                      path_C8_info = path_C8_info, 
                                                      path_C18_info = path_C18_info, 
                                                      path_HILIC_info = path_HILIC_info, 
                                                      QC_label = "QC-pooled_ref", 
                                                      Metabs_long = Metabs_long)),

#save CV file for use
tar_target(Metabolite_CV_and_missingness, make_metabolite_QC$QC_info),
tar_target(Metabolite_CV_and_missingness_filename, paste0("MESA_TOPMed_Metabolite_QCfile_", Sys.Date(), ".csv")),


#export CV as csv
tar_target(Metabolite_CV_and_missingness_csv,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, Metabolite_CV_and_missingness_filename)
             readr::write_csv(Metabolite_CV_and_missingness, out_path)
             out_path
           },
           format = "file"
),

#----------------------------------------------------------------------------#
#--------------------------------4. Duplicates flag file---------------------#
#----------------------------------------------------------------------------#

tar_target(Duplicates_flag_file, flag_duplicates_function(QC_file = Metabolite_CV_and_missingness, 
                                                          mapping_file = Mapping_file)),

#Save
#save duplicate flagging file for use
tar_target(Duplicates_flag_file_filename, paste0("MESA_TOPMed_Metabolite_Duplicateflag_", Sys.Date(), ".csv")),


#export CV as csv
tar_target(Duplicates_flag_csv,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, Duplicates_flag_file_filename)
             readr::write_csv(Duplicates_flag_file, out_path)
             out_path
           },
           format = "file"
),
           
#----------------------------------------------------------------------------#
#--------------------------------4. Make_final clean file---------------------#
#----------------------------------------------------------------------------#

tar_target(Final_metabs_long_outfile, make_final_metabs_final_function(duplicate_flagging_file = Duplicates_flag_file, 
                                                               metabolite_file = Metabs_long)),
tar_target(Final_metabs_long, Final_metabs_long_outfile$metabolite_file_nodupes),
tar_target(Final_metabs_long_info, Final_metabs_long_outfile$metabolite_file_nodupes_dims),

#Save
#save duplicate flagging file for use
tar_target(Final_metabs_long_filename, paste0("MESA_TOPMed_Metabolite_cleanfile_", Sys.Date(), ".csv")),


#export CV as csv
tar_target(Final_metabs_long_csv,
           {
             out_dir  <- "outputs"
             dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
             out_path <- file.path(out_dir, Final_metabs_long_filename)
             readr::write_csv(Final_metabs_long, out_path)
             out_path
           },
           format = "file"
),


#----------------------------------------------------------------------------#
#--------------------------------Quarto file---------------------------------#
#----------------------------------------------------------------------------#

tar_quarto(
  build_proteins_quarto,
  path = "/media/Analyses/Make-MESA-TOPMed-Metabolite-Table/Make-MESA-TOPMed-Metabolite-Table-README.qmd",
  quiet = FALSE
)

)