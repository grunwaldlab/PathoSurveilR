#' Find the first directory containing pathogensurveillance output
#'
#' Recursively searches a directory to find the first directory (starting with
#' the input directory) that contains a file called
#' '.pathogensurveillance_output.yml'.
#'
#' @param path `character`. The starting paths to search. Defaults to current
#'   directory.
#' @param allow_nested `logical`. If `FALSE`, don't search subdirectories of
#'   directories that already contain the target file. If `TRUE`, continue
#'   searching all subdirectories.
#'
#' @return `character` vector of paths to directories containing the target
#'   file.
#'
#' @examples
#' \dontrun{
#' # Search current directory and subdirectories (non-nested)
#' find_output_dirs()
#'
#' # Search a specific path with nested searching allowed
#' find_output_dirs("~/projects", allow_nested = TRUE)
#' }
#'
#' @export
find_output_dirs <- function(path = ".", allow_nested = FALSE) {
  # Normalize the path
  path <- normalizePath(path, mustWork = TRUE)
  
  # Internal recursive function that accumulates results
  recursive_search <- function(current_path) {
    # Check if target file exists in current directory
    if (is_output_directory(current_path)) {
      found <- current_path
      
      # If we're not allowing nested searches, return immediately
      if (!allow_nested) {
        return(found)
      }
    } else {
      found <- character(0)
    }
    
    # Search each subdirectory
    dirs <- list.dirs(current_path, full.names = TRUE, recursive = FALSE)
    for (dir in dirs) {
      subdir_results <- .search_dir(dir)
      found <- c(found, subdir_results)
    }
    
    return(found)
  }
  
  # Recursively search all input paths and combine the results
  unlist(lapply(path, recursive_search))
}


#' @keywords internal
is_output_directory <- function(path) {
  pipeline_yml_path <- file.path(path, '.pathogensurveillance_output.json')
  file.exists(pipeline_yml_path)
}


#' Find paths for specific pathogensurveillance output files
#'
#' Find paths for specific pathogensurveillance output files and associated
#' metadata
#'
#' @param outdir_path The path to an output directory of pathogensurviellance
#' @param targets The names of one or more output types to search for.
#' @param simplify If `TRUE`, combine all of the output found into a single
#'   table using the `long` option to decide how to format it.
#' @param long If `TRUE` when combining output into a single table, put all
#'   paths in a single `path` column rather than having a column for each
#'   `target`.
#'
#' @export
find_path_data <- function(outdir_path, targets = NULL, simplify = TRUE, long = TRUE) {
  # Check that outdir_path is valid pathogensurveillance output
  if (!is_output_directory(outdir_path)) {
    stop('`outdir_path` does not seem to a pathogensurveillance output directory.')
  }
  outdir_path <- normalizePath(outdir_path)
  
  # Load output schema metadata
  metadata <- parse_output_meta_json(outdir_path)
  
  # Return data on all available data if targets is undefined
  if (is.null(targets)) {
    targets <- names(metadata$outputs)
  }
  
  # Search for files
  path_data <- lapply(targets, function(target) {
    target_meta <- metadata$outputs[[target]]
    do.call(rbind, lapply(target_meta$sources, function(source) {
      source_path <- file.path(outdir_path, source$path)
      possible_paths <- list.files(path = source_path, recursive = TRUE)
      matching_paths <- possible_paths[grepl(possible_paths, pattern = source$pattern)]
      matching_paths <- matching_paths[file.exists(file.path(source_path, matching_paths))]
      if (length(matching_paths) == 0) {
        return(NULL)
      }
      capture_group_key <- rep(list(character(0)), length(source$capture_groups))
      names(capture_group_key) <- source$capture_groups
      output <- utils::strcapture(matching_paths, pattern = source$pattern, proto = capture_group_key)
      output[names(source$constants)] <- source$constants
      output$path <- file.path(source_path, matching_paths)
      tibble::as_tibble(output)
    }))
  })
  names(path_data) <- targets
  
  # Combine the results for each output type into a single table
  if (simplify) {
    path_data <- path_data[! vapply(path_data, FUN.VALUE = logical(1), is.null)]
    extra_columns <- unique(unlist(lapply(path_data, colnames)))
    extra_columns <- extra_columns[! extra_columns %in% targets]
    extra_columns <- extra_columns[extra_columns != 'path']
    path_data <- lapply(path_data, function(x) {
      unused_cols <- extra_columns[! extra_columns %in% colnames(x)]
      x[unused_cols] <- rep(list(NA_character_), length(unused_cols))
      x
    })
    if (long) {
      path_data <- lapply(names(path_data), function(n) {
        x <- path_data[[n]]
        x$target <- n
        x
      })
      path_data <- do.call(combine_data_frames, path_data)
    } else {
      path_data <- lapply(names(path_data), function(n) {
        x = path_data[[n]]
        colnames(x)[colnames(x) == 'path'] <- n
        x
      })
      path_data <- Reduce(path_data, f = function(x, y) {
        merge(x, y, by = extra_columns, all = TRUE)
      })
      path_data <- tibble::as_tibble(path_data)
    }
  }

  return(path_data)
}


#' @keywords internal
parse_output_meta_json <- function(outdir_path) {
  # Check that outdir_path is valid pathogensurveillance output
  if (!is_output_directory(outdir_path)) {
    stop('`outdir_path` does not seem to a pathogensurveillance output directory.')
  }
  outdir_path <- normalizePath(outdir_path)
  
  meta_path <- file.path(outdir_path, '.pathogensurveillance_output.json')
  RcppSimdJson::fload(meta_path, max_simplify_lvl = 'vector')
}

