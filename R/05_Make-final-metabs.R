make_final_metabs_final_function <- function(duplicate_flagging_file, metabolite_file)
{
  
  keep_metabs <- subset(duplicate_flagging_file, Retain== "Unique or missing HMDB ID" | Retain== "Duplicated HMDB ID with lowest CV")
    metabolite_file_nodupes <- metabolite_file[,names(metabolite_file)=="sidno" | 
                                                 names(metabolite_file)=="idno" |
                                                 names(metabolite_file)=="TOM_ID" |
                                                 names(metabolite_file)=="subject_id" |
                                                 names(metabolite_file)=="exam" | names(metabolite_file) %in% keep_metabs$Metabolite]
  

#----------------------------Output----------------------------------------#
  
   list( 
     metabolite_file_nodupes = metabolite_file_nodupes,
     metabolite_file_nodupes_dims = dim(metabolite_file_nodupes)
   )
}