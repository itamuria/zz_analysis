# Rename SCE clusters to population
# From excel file. Structure: Poblacion, cluster
# Example:
# uncorrected <- renames_clusters_deg (sce_unc = uncorrected,
                # unc_cl_name = "tcell_clustering",
                # new_var = "mainPop", df_table = df_table,
                # var_name = "Poblacion", cluster_name = "cluster")

renames_clusters_deg <- function(sce_unc,
                                 unc_cl_name,
                                 new_var, df_table,
                                 var_name, cluster_name)
{
  
  tab_names <- df_table
  
  # how many populations
  pop_names <- names(table(tab_names[,var_name]))
  
  colData(sce_unc)[new_var] <- 0
  
  new_vec <- rep("vacio",nrow(colData(sce_unc)))
  
  for(p in pop_names)
  {
    print(p)
    temp <- tab_names[tab_names[,var_name] %in% p, ]
    
    new_vec <- ifelse(unlist(colData(sce_unc)[unc_cl_name] %in% temp[,cluster_name]), p, new_vec)
    
  }
  
  colData(sce_unc)[new_var] <- new_vec
  
  return(sce_unc)
  
}