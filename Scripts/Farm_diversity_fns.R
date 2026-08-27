# Shared helpers for the Farm_diversity pipeline ----

## Return the most recently modified file in `dir` whose name matches `pattern`; lets the matching / modelling scripts pick up the latest date-stamped Tax_div_* export from 02_Analysis_iNEXT.qmd without hard-coding the date
latest_file <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) {
    stop("No file matching '", pattern, "' in ", dir, call. = FALSE)
  }
  hits[which.max(file.mtime(hits))]
}
