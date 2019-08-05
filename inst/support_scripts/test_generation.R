cat(paste0("test_that(\"", 
           apply(
             expand.grid(c("single", "multiple"), 
                         c("short", "medium-length", "long"), 
                         c("ambiguous", "non-ambiguous"), 
                         c("nucleotides", "aminoacids"), 
                         c("", "with specified type", "with specified is_clean", "with specified both type and is_clean"),
                         stringsAsFactors = FALSE), 1, function(row) 
                           paste(row, sep = " ", collapse = " ")), 
           " sequences constructing\", {})\n"))

cat(paste0("test_that(\"", 
           apply(
             expand.grid(c("ambiguous", "non-ambiguous"), 
                         c("nucleotides", "aminoacids"), 
                         c("with specified type", "with specified is_clean", "with specified both type and is_clean"),
                         stringsAsFactors = FALSE), 1, function(row) 
                           paste(row, sep = " ", collapse = " ")), 
           " sequences constructing\", {})\n"))
