# Function to clean and update marker names for a single cell type list
update_marker_names <- function(marker_list, seurat_obj) {
  
  # Use map to iterate over each cell type vector in the list
  updated_list <- map(marker_list, function(marker_vector) {
    
    # Use map_chr to iterate over each individual marker in the vector
    map_chr(marker_vector, function(marker) {
      
      # 1. Find gene name (case-insensitive)
      new_name <- grep(marker, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
      
      # 2. Check and return the appropriate name
      if (length(new_name) == 1) {
        return(new_name) # Found exactly one match, use the official name
      } else if (length(new_name) > 1) {
        # Multiple hits found. The original code's logic here was flawed (if (marker %in% marker)),
        # but the intent seems to be: if we can't be sure, keep the original name to avoid errors.
        # We'll print a warning and keep the original name.
        # warning(glue("Multiple hits found for '{marker}'. Keeping original name."))
        return(marker) 
      } else {
        # No matches found. Print an error message and keep the original name.
        # Keeping the original name allows you to manually fix it later.
        print(glue("'{marker}' not found in Seurat object."))
        return(marker) 
      }
    })
  })
  
  # Assign the corrected list back to the original names
  names(updated_list) <- names(marker_list)
  return(updated_list)
}