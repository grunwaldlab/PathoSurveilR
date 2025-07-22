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
#' @param outdir_path One or more paths to output directories of specific output
#'   files of pathogensurviellance
#' @param target The names of one or more output types to search for. Setting to
#'   `NULL` causes data for all available paths to be returned.
#' @param simplify If `TRUE`, combine all of the output found into a single
#'   table using the `long` option to decide how to format it.
#' @param long If `TRUE`, put all paths in a single column rather than having
#'   each target in their own column.
#' @param all_metadata Include columns for all metadata regarding the discovery
#'   and parsing of files.
#' @param must_exist If `TRUE`, an error is generated if the specified target
#'   does not exist. Otherwise `NULL` is returned for missing `target`s
#'
#' @export
find_ps_paths <- function(path, target, simplify = TRUE, long = TRUE, all_metadata = FALSE, must_exist = ! is.null(target)) {
  
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
      my_target <- names(metadata$outputs)
    } else {
      my_target <- target
    }
    
    # Search for files
    path_data <- lapply(my_target, function(t) {
      target_meta <- metadata$outputs[[t]]
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
    names(path_data) <- my_target
    
    # Filter possible paths by original path supplied
    path_data <- lapply(path_data, function(x) {
      if (is.null(x)) {
        return(NULL)
      } else {
        out <- x[startsWith(x$path, one_path), , drop = FALSE]
        if (nrow(out) == 0) {
          return(NULL)
        } else {
          return(out)
        }
      }
    })
    
    return(path_data)
  })
  all_path_data <- unlist(all_path_data, recursive = FALSE)
  
  # If `target` is specified, error if any targets are not found
  if (must_exist) {
    is_missing <- vapply(all_path_data, FUN.VALUE = logical(1), is.null)
    missing_targets <- names(all_path_data)[is_missing]
    if (length(missing_targets) > 0) {
      valid_targets <- known_ps_outputs(outdir_path = path, exists = TRUE)$target
      suggestions <- unlist(lapply(missing_targets, suggest_option, options = valid_targets))
      stop(
        call. = FALSE,
        'The following targets could not be found:\n  ', paste0('"', missing_targets, '"', collapse = ', '), '\n',
        ifelse(length(suggestions) == 0, '', paste0('Valid targets with similar names include:\n  ', paste0('"', suggestions, '"', collapse = ', '), '\n')),
        'To see all valid targets for this output directory use the following:\n',
        '  print_ps_outputs("', path, '", exists = TRUE)\n' 
      )
    }
  }
  

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
  unique_paths <- unique_paths[order(nchar(unique_paths))]
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
  lapply(strsplit(path, split = '/'), function(split_path) {
    vapply(seq_len(length(split_path)), FUN.VALUE = character(1), function(n) {
      paste0(split_path[1:n], collapse = '/')
    })
  })
  # lapply(path, function(p) {
  #   output <- p
  #   while (! p %in% c('/', '.', '..')) {
  #     p <- dirname(p)
  #     output <- c(output, p)
  #   }
  #   return(output)
  # })
}


#' Prints available pipeline output files
#'
#' Prints available pipeline output files grouped by type
#'
#' @inheritParams known_ps_outputs
#'
#' @return Invisibly returns a nested list of categories with their items and descriptions.
#'   Primarily called for its side effect of printing formatted output to the console.
#'
#' @export
print_ps_outputs <- function(outdir_path, exists = TRUE) {
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
known_ps_outputs <- function(outdir_path, exists = TRUE) {
  # Load output schema metadata
  metadata <- parse_output_meta_json(outdir_path)
  
  # Filter by output type that have file associated with them
  if (exists) {
    all_path_data <- find_ps_paths(outdir_path, target = NULL, simplify = TRUE, long = TRUE)
    all_output_types <- unique(all_path_data$target)
    metadata$outputs <- metadata$outputs[names(metadata$outputs) %in% all_output_types]
  }
  
  desc_data <- do.call(rbind, lapply(names(metadata$outputs), function(output_name) {
    output <- metadata$outputs[[output_name]]
    description <- output$description
    category <- output$category
    parser <- output$parsers$r
    combiner <- output$combiners$r
    if (is.null(parser)) {
      parser <- NA_character_
    }
    if (is.null(combiner)) {
      combiner <- NA_character_
    }
    if (is.null(category) || category == "") {
      category <- 'Miscellaneous'
    }
    if (is.null(description)) {
      category <- ''
    }
    
    # Process each source in the output
    sources_df <- do.call(rbind, lapply(output$sources, function(source) {
      path <- source$path
      pattern <- source$pattern
      
      # Create name_schema by replacing capture groups
      if (!is.null(source$capture_groups) && length(source$capture_groups) > 0) {
        # Split the pattern at capture groups
        pattern_parts <- strsplit(pattern, "\\(.+?\\)", perl = TRUE)[[1]]
        capture_groups <- paste0("<", source$capture_groups, ">")
        
        # Interleave pattern parts and capture groups
        # (This works when capture groups are at the end)
        name_schema <- paste0(
          paste(
            c(rbind(pattern_parts, c(capture_groups, ""))), 
            collapse = ""
          )
        )
        
        # Clean up any regex characters
        name_schema <- gsub("\\^|\\$|\\\\|\\?", "", name_schema)
      } else {
        name_schema <- gsub("\\^|\\$|\\\\|\\?", "", pattern)
      }
      
      # Move parts of the file name schema that are directories into the path
      split_schema <- strsplit(name_schema, split = '/')[[1]]
      name_schema <- split_schema[[length(split_schema)]]
      path <- paste0(c(path, split_schema[seq_len(length(split_schema) - 1)]), collapse = '/')
      
      # Combine path and name_schema
      schema <- paste0(path, '/', name_schema)
      
      data.frame(
        target = output_name,
        description = description,
        category = category,
        parser = parser,
        combiner = combiner,
        schema = schema,
        stringsAsFactors = FALSE
      )
    }))
    
    sources_df
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


#' Print pathogensurveillance outdir schema
#'
#' Print the output directory layout for a given pipeline output
#'
#' @inheritParams known_ps_outputs
#'
#' @export
print_outdir_schema <- function(outdir_path, exists = TRUE) {
  
  path_data <- known_ps_outputs(outdir_path, exists = exists)
  
  # Sort alphabetically by path
  path_data <- path_data[order(path_data$schema), , drop = FALSE]
  
  # Split paths into components
  paths <- strsplit(path_data$schema, "/", fixed = TRUE)
  
  # Create a data structure to track directory hierarchy
  dir_levels <- list()
  max_depth <- max(sapply(paths, length))
  
  # Initialize the tree structure
  tree <- data.frame(
    name = ".",
    depth = 0,
    is_last = FALSE,
    parent = NA,
    description = NA,
    stringsAsFactors = FALSE
  )
  
  # Process each path
  for (i in seq_along(paths)) {
    current_parent <- 1  # Start from root
    for (depth in seq_along(paths[[i]])) {
      component <- paths[[i]][depth]
      is_file <- (depth == length(paths[[i]]))
      
      # Check if this component already exists at this level
      existing <- which(
        tree$parent == current_parent & 
          tree$name == component & 
          tree$depth == depth
      )
      
      if (length(existing) == 0) {
        # Add new node
        new_node <- data.frame(
          name = component,
          depth = depth,
          is_last = FALSE,  # Will update later
          parent = current_parent,
          description = ifelse(is_file, path_data$description[i], NA),
          stringsAsFactors = FALSE
        )
        tree <- rbind(tree, new_node)
        current_parent <- nrow(tree)
      } else {
        current_parent <- existing[1]
      }
    }
  }
  
  # Determine which nodes are last in their groups
  for (parent_id in unique(tree$parent)) {
    if (!is.na(parent_id)) {
      children <- which(tree$parent == parent_id)
      tree$is_last[children[length(children)]] <- TRUE
    }
  }
  
  # Prepare indentation patterns and calculate line lengths
  indent_symbols <- character(nrow(tree))
  line_lengths <- integer(nrow(tree))
  max_line_length <- 0
  
  for (i in 2:nrow(tree)) {
    if (tree$depth[i] == 1) {
      indent_symbols[i] <- ifelse(tree$is_last[i], "└── ", "├── ")
    } else {
      # Get all ancestors
      ancestors <- c()
      current <- i
      while (!is.na(tree$parent[current])) {
        ancestors <- c(tree$parent[current], ancestors)
        current <- tree$parent[current]
      }
      
      # Build indentation
      indent <- ""
      for (a in ancestors[-1]) {  # Skip root
        if (tree$is_last[a]) {
          indent <- paste0(indent, "    ")
        } else {
          indent <- paste0(indent, "│   ")
        }
      }
      indent_symbols[i] <- paste0(
        indent,
        ifelse(tree$is_last[i], "└── ", "├── ")
      )
    }
    
    # Calculate line length (without description)
    line_lengths[i] <- nchar(indent_symbols[i]) + nchar(tree$name[i])
    if (!is.na(tree$description[i])) {
      # Update max line length if this is a file with description
      max_line_length <- max(max_line_length, line_lengths[i])
    }
  }
  
  # Print the tree
  cat(".\n")
  for (i in 2:nrow(tree)) {  # Skip root
    line <- paste0(indent_symbols[i], tree$name[i])
    if (!is.na(tree$description[i])) {
      # Calculate padding needed to align descriptions
      padding <- max_line_length - line_lengths[i] + 2
      line <- paste0(line, strrep(" ", padding), ": ", tree$description[i])
    }
    cat(line, "\n")
  }
}
