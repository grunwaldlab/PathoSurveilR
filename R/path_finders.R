#' Find directories containing pathogensurveillance output
#'
#' Recursively searches a directory to find the first directory (starting with
#' the input directory) that contains a file called
#' '.pathogensurveillance_output.json'.
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
find_outdirs <- function(path = ".", allow_nested = FALSE) {
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
#' @param target The names of one or more output types to search for.
#' @param simplify If `TRUE`, combine all of the output found into a single
#'   table using the `long` option to decide how to format it.
#' @param long If `TRUE`, put all paths in a single column rather than having
#'   each target in their own column.
#' @param all_metadata Include columns for all metadata regarding the discovery
#'   and parsing of files.
#'
#' @export
find_ps_paths <- function(path, target = NULL, simplify = TRUE, long = TRUE, all_metadata = FALSE) {
  
  all_path_data <- lapply(path, function(one_path) {
    # Search for the first parent path that is a pathogensurveillance output directory
    outdir_path <- get_ps_outdir_from_path(one_path)
    if (is.na(outdir_path)) {
      stop('`path` value "', one_path, '" does not seem to be within a pathogensurveillance output directory.')
    }
    
    # Save the leftover part of the path to filter the results with later
    path_suffix <- substr(one_path, start = nchar(outdir_path) + 2, stop = nchar(one_path))
    
    # Load output schema metadata
    metadata <- parse_output_meta_json(outdir_path)
    
    # Return data on all available data if target is undefined
    if (is.null(target)) {
      target <- names(metadata$outputs)
    }
    
    # Search for files
    path_data <- lapply(target, function(target) {
      target_meta <- metadata$outputs[[target]]
      do.call(rbind, lapply(target_meta$sources, function(source) {
        source_path <- file.path(outdir_path, source$path)
        possible_paths <- list.files(path = source_path, recursive = TRUE,
                                     all.files = TRUE, include.dirs = TRUE)
        matching_paths <- possible_paths[grepl(possible_paths, pattern = source$pattern)]
        matching_paths <- matching_paths[file.exists(file.path(source_path, matching_paths))]
        if (length(matching_paths) == 0) {
          return(NULL)
        }
        
        # Parse columns derived from values embedded in the file path
        if (is.null(source$capture_groups)) {
          output <- tibble::tibble(path = file.path(source_path, matching_paths))
        } else {
          capture_group_key <- rep(list(character(0)), length(source$capture_groups))
          names(capture_group_key) <- source$capture_groups
          output <- utils::strcapture(matching_paths, pattern = source$pattern, proto = capture_group_key)
          output$path <- file.path(source_path, matching_paths)
        }
        
        # Add input parsing metadata
        output$outdir_path <- outdir_path
        if (all_metadata) {
          if (is.null(target_meta$parsers$r)) {
            output$parser <- NA_character_
          } else {
            output$parser <- target_meta$parsers$r
          }
          if (is.null(target_meta$combiners$r)) {
            output$combiner <- NA_character_
          } else {
            output$combiner <- target_meta$combiners$r
          }
          output$category <- target_meta$category
          output$description <- target_meta$description
        }
        
        # Add constants values
        output[names(source$constants)] <- source$constants
        
        tibble::as_tibble(output)
      }))
    })
    names(path_data) <- target
    path_data <- path_data[! sapply(path_data, is.null)]
    
    # Filter possible paths by original path supplied
    path_data <- lapply(path_data, function(x) {
      x[startsWith(x$path, one_path), , drop = FALSE]
    })
    path_data <- path_data[sapply(path_data, nrow) > 0]
    
    return(path_data)
  })
  all_path_data <- unlist(all_path_data, recursive = FALSE)

  # Simplify output to a single table
  if (simplify) {
    # Remove empty results
    all_path_data <- all_path_data[! vapply(all_path_data, FUN.VALUE = logical(1), is.null)]
    
    # Add placeholder columns so that all tables have the same columns
    extra_columns <- unique(unlist(lapply(all_path_data, colnames)))
    extra_columns <- extra_columns[! extra_columns %in% names(all_path_data)]
    extra_columns <- extra_columns[! extra_columns %in% c('path', 'category', 'parser', 'description')]
    all_path_data <- lapply(all_path_data, function(x) {
      unused_cols <- extra_columns[! extra_columns %in% colnames(x)]
      x[unused_cols] <- rep(list(NA_character_), length(unused_cols))
      x
    })
    
    # Combine into a single table
    if (long) {
      all_path_data <- lapply(seq_len(length(all_path_data)), function(i) {
        tibble::as_tibble(cbind(target = names(all_path_data)[i], all_path_data[[i]]))
      })
      all_path_data <- combine_data_frames(all_path_data)
    } else {
      all_path_data <- lapply(seq_len(length(all_path_data)), function(i) {
        x = all_path_data[[i]]
        colnames(x)[colnames(x) == 'path'] <- names(all_path_data)[i]
        x
      })
      all_path_data <- Reduce(all_path_data, f = function(x, y) {
        merge(x, y, by = extra_columns, all = TRUE)
      })
      all_path_data <- tibble::as_tibble(all_path_data)
    }
  }  

  return(all_path_data)
}


#' Get parent directory that is a pathogensurveillance output directory
#' 
#' @param path A path to a file/folder within a pathogensurveillance output directory
#'
#' @export
get_ps_outdir_from_path <- function(path) {
  # Check all input paths and parent directories
  path_parts <- split_path(path)
  unique_paths <- unique(unlist(path_parts))
  is_out_key <- is_output_directory(unique_paths)
  names(is_out_key) <- unique_paths
  
  # Return most specific output directory for each path
  vapply(path_parts, FUN.VALUE = character(1), function(x) {
    if (any(is_out_key[x])) {
      return(x[max(which(is_out_key[x]))])
    } else {
      return(NA_character_)
    }
  })
}


#' @keywords internal
split_path <- function(path) {
  lapply(path, function(p) {
    output <- p
    while (! p %in% c('/', '.', '..')) {
      p <- dirname(p)
      output <- c(output, p)
    }
    return(output)
  })
}


#' Prints available pipeline output files
#'
#' Prints available pipeline output files grouped by type
#'
#' @inheritParams known_outputs
#'
#' @return Invisibly returns a nested list of categories with their items and descriptions.
#'   Primarily called for its side effect of printing formatted output to the console.
#'
#' @export
print_ps_outputs <- function(outdir_path, exists = FALSE) {
  desc_data <- known_ps_outputs(outdir_path, exists = exists)

  # Get print text for each line
  desc_data$print_text <- paste0(
    '  ',
    format(desc_data$target, width = max(nchar(desc_data$target)), justify = "left"),
    ' : ',
    desc_data$description
  )
  
  # Print results grouped by category
  category_print_text <- lapply(split(desc_data, desc_data$category), function(category_data) {
    category <- unique(category_data$category)
    width <- max(nchar(desc_data$print_text))
    header <- paste0(
      format(category, width = width, justify = 'left'), '\n',
      format(strrep("\U00AF", width), width = width, justify = 'left')
    )
    paste0(c(header, category_data$print_text), collapse = '\n')
  })
  cat(paste0(category_print_text, collapse = '\n\n'))
}


#' Get which outputs are available
#'
#' Get which outputs are available for a given pathogensurveillance output
#' directory.
#'
#' @param outdir_path The path to an output directory of pathogensurviellance
#' @param exists If `TRUE`, information is only returned for output type that
#'   exist, rather than all output type that are known.
#'
#' @export
known_ps_outputs <- function(outdir_path, exists = FALSE) {
  # Load output schema metadata
  metadata <- parse_output_meta_json(outdir_path)
  
  # Filter by output type that have file associated with them
  if (exists) {
    all_path_data <- find_path_data(outdir_path, simplify = TRUE, long = TRUE)
    all_output_types <- unique(all_path_data$target)
    metadata$outputs <- metadata$outputs[names(metadata$outputs) %in% all_output_types]
  }
  
  # Organize outputs by category
  desc_data <- do.call(rbind, lapply(names(metadata$outputs), function(key) {
    item <- metadata$outputs[[key]]
    category <- item$category
    description <- item$description
    parser <- item$parser$r
    if (is.null(parser)) {
      parser <- NA_character_
    }
    if (is.null(category) || category == "") {
      category <- 'Miscellaneous'
    }
    if (is.null(description)) {
      category <- ''
    }
    data.frame(
      target = key,
      description = description,
      category = category,
      parser = parser,
      stringsAsFactors = FALSE
    )
  }))
  
  tibble::as_tibble(desc_data)
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

