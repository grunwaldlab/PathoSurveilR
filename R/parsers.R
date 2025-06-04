#' Get summary of pipeline status data
#'
#' Return a [tibble::tibble()] (table) with the numbers of issues encountered by the
#' pipeline. The contents of all status message files found in the given paths
#' will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with counts of messages attributes if
#'   `simplify = TRUE` or a list of such [tibble::tibble()]s for each output
#'   directory found if `simplify = FALSE`.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_parsed_summary(path)
#'
#' @export
status_message_summary <- function(path) {
  message_data <- find_ps_data(path, target = 'messages', simplify = TRUE)
  if (nrow(message_data) == 0) {
    output <- tibble::tibble(
      workflow = character(0),
      errors =  numeric(0),
      warnings = numeric(0),
      notes = numeric(0),
      samples = numeric(0),
      references = numeric(0),
      report_groups = numeric(0)
    )
  } else {
    output <- do.call(rbind, lapply(split(message_data, message_data$workflow), function(subset) {
      tibble::tibble(
        workflow = unique(subset$workflow),
        errors = sum(subset$level == 'ERROR'),
        warnings = sum(subset$level == 'WARNING'),
        notes = sum(subset$level == 'NOTE'),
        samples = length(unique(subset$sample_id)),
        references = length(unique(subset$reference_id)),
        report_groups = length(unique(subset$report_id))
      )
    }))
  }
  return(output)
}


#' Get parsed sendsketch results
#'
#' Return a [tibble::tibble()] (table) with the results from sendsketch. The
#' contents of all sendsketch files found in the given paths will be combined.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param sample_id One or more sample IDs to parse the results for. Default:
#'   include all samples.
#' @param only_best Only return the best hit for each combination of report
#'   group and sample.
#' @param only_shared If `TRUE`, only return the ranks that are present in all
#'   of the inputs. Only has an effect if `update_taxonomy` is `TRUE`.
#' @param update_taxonomy If `TRUE` look up the current taxonomy classification
#'   using the NCBI taxon ID rather than using the existing one and add columns for each rank. 
#' @inheritParams postprocess_table_list
#'
#' @return A [tibble::tibble()] with the sendsketch output combined
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_parsed(path)
#'
#' @export
sendsketch_parsed <- function(path, sample_id = NULL, only_best = FALSE, only_shared = FALSE, update_taxonomy = FALSE, simplify = TRUE) {
  path_data <- find_ps_paths(path, target = 'sendsketch', simplify = FALSE)[[1]]
  
  if (! is.null(sample_id)) {
    path_data <- path_data[path_data$sample_id %in% sample_id]    
  }
  
  # Parse all files and combine them into one data frame
  if (is.null(path_data) || nrow(path_data) == 0) {
    sketch_data <- make_empty_data_frame(
      c("sample_id", "WKID", "KID", "ANI", "SSU", "SSULen", "Complt", 
        "Contam", "Contam2", "uContam", "Score", "E-Val", "Depth", "Depth2", 
        "Volume", "RefHits", "Matches", "Unique", "Unique2", "Unique3", 
        "noHit", "Length", "TaxID", "ImgID", "gBases", "gKmers", "gSize", 
        "gSeqs", "GC", "rDiv", "qDiv", "rSize", "qSize", "cHits", "taxName", 
        "file", "seqName", "taxonomy", "outdir_path")
    )
  } else {
    raw_data <- find_ps_data(path, target = 'sendsketch') # TODO: change out for function that only parses the required datasets
    sketch_data <- do.call(rbind, lapply(seq_len(nrow(path_data)), function(i) {
      data <- raw_data[[path_data$path[i]]]
      return(cbind(
        sample_id = rep(path_data$sample_id[i], nrow(data)),
        data
      ))
    }))
    sketch_data <- tibble::as_tibble(sketch_data)
  }
  
  # Convert columns to the correct type
  sketch_data[] <- lapply(sketch_data, as.character)
  numeric_cols <- c(6, 11:22, 28:34)
  sketch_data[numeric_cols] <- lapply(sketch_data[numeric_cols], as.numeric)
  sketch_data$WKID <- as.numeric(gsub("%", "", sketch_data$WKID))
  sketch_data$ANI <- as.numeric(gsub("%", "", sketch_data$ANI))
  sketch_data$Complt <- as.numeric(gsub("%", "", sketch_data$Complt))
  
  # Filter for top hits
  if (only_best && nrow(sketch_data) > 0) {
    sketch_data <- sendsketch_best_hits(sketch_data)
  }

  # Look up taxonomy based on taxon ID for most up to date taxonomy
  if (update_taxonomy && nrow(sketch_data) > 0)  {
    tax_data <- sendsketch_taxonomy_parsed(input = sketch_data, only_shared = only_shared, simplify = TRUE)
    
    # Combine taxonomy data based on taxon ID while preserving column order
    original_columns <- colnames(sketch_data)
    sketch_data <- merge(x = sketch_data, y = tax_data, by.x = 'TaxID', by.y = 'taxon_id', all.x = TRUE)
    new_columns <- colnames(sketch_data)[! colnames(sketch_data) %in% original_columns]
    column_order <- c(original_columns, new_columns)
    sketch_data <- sketch_data[column_order]
  }
  
  tibble::as_tibble(sketch_data)
}


#' Lookup classifications using taxon IDs
#' 
#' @param match_input If `TRUE`, reduplicate the output to match 1:1 with the input.
#' @param simplify If `TRUE`, combine list of tibbles into a single tibble
#' 
#' @keywords internal
lookup_ncbi_taxon_id_classification <- function(taxon_ids, match_input = TRUE, simplify = FALSE) {
  all_tax_ids <- unique(taxon_ids)
  tax_id_batches <- split(all_tax_ids, ceiling(seq_along(all_tax_ids) / 50))
  output <- lapply(tax_id_batches, function(id_batch) {
    # Look up taxon info for each ID
    raw_xml = rentrez::entrez_fetch(db = 'taxonomy', id = id_batch, rettype = 'xml')
    Sys.sleep(0.33)
    
    # Parse XML to get IDs, names, and ranks of parent taxa
    xml_data <- xml2::read_xml(raw_xml)
    result_data <- xml2::xml_find_all(xml_data, "/TaxaSet/Taxon")
    get_one_value <- function(key) {
      lapply(result_data, function(x) {
        query_value <- xml2::xml_text(xml2::xml_find_all(x, paste0("./", key)))
        lineage_values <- xml2::xml_text(xml2::xml_find_all(x, paste0("./LineageEx/Taxon/", key)))
        c(lineage_values, query_value)
      })
    }
    taxon_ids <- get_one_value('TaxId')
    taxon_names <- get_one_value('ScientificName')
    taxon_ranks <- get_one_value('Rank')
    
    # Format as table
    lapply(seq_along(id_batch), function(i) {
      tibble::tibble(
        query_id = id_batch[i],
        id = taxon_ids[[i]],
        name = taxon_names[[i]],
        rank = taxon_ranks[[i]]
      )
    })
  })
  output <- unlist(output, recursive = FALSE)
  names(output) <- all_tax_ids
  
  if (match_input) {
    output <- output[taxon_ids]
  }
  
  if (simplify) {
    output <- do.call(rbind, output)
  }
  
  return(output)
}


#' Get parsed sendsketch taxonomy
#'
#' Return the taxonomic classification in sendsketch results associated with
#' directory paths containing `pathogensurveillance` output.
#'
#' @param input The path to one or more folders that contain pathogensurveillance
#'   output or a table in the form of [sendsketch_parsed()] output.
#' @param ranks_as_cols If `TRUE`, return a table with a column for each rank
#'   instead of a table with one taxon per row.
#' @param include_sample_id description
#' @param only_shared If `TRUE`, only return the ranks that are present in all
#'   of the inputs. Only has an effect if `as_table` is `TRUE`.
#' @inheritParams sendsketch_parsed
#' @inheritParams postprocess_table_list
#'
#' @return A [base::character()] vector of taxonomic classifications, each
#'   delimited with `;`, named by sample IDs.
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_taxonomy_parsed(path)
#'
#' @export
sendsketch_taxonomy_parsed <- function(input, ranks_as_cols = TRUE, only_shared = FALSE,
                                       sample_id = NULL, only_best = FALSE, simplify = TRUE) {
  
  # Used existing data if supplied as a table, otherwise find and parse data based on a path
  if (inherits(input, 'data.frame')) {
    sendsketch_data <- input
  } else {
    sendsketch_data <- sendsketch_parsed(input, sample_id = sample_id, only_best = only_best,
                                         update_taxonomy = FALSE, simplify = FALSE)
  }

  class_data <- lookup_ncbi_taxon_id_classification(sendsketch_data$TaxID, match_input = FALSE, simplify = FALSE)
  
  if (nrow(sendsketch_data) == 0) {
    if (as_table) {
      return(tibble::as_tibble(make_empty_data_frame(
        c("taxon_id", "phylum", "order", "family", "genus", "species")
      )))
    } else {
      return(tibble::as_tibble(make_empty_data_frame(
        c("query_id", "id", "name", "rank")
      )))
    }
  }
  
  unique_tax_ids <- unique(sendsketch_data$TaxID)
  table_class_data <- class_data[unique_tax_ids]
  
  # Remove data for ranks not shared among all results
  if (only_shared) {
    all_ranks <- unlist(lapply(table_class_data, function(x) unique(x$rank)))
    is_shared <- table(all_ranks) == length(table_class_data)
    shared_ranks <- names(is_shared)[is_shared]
    table_class_data <- lapply(table_class_data, function(x) x[x$rank %in% shared_ranks, , drop = FALSE])
  }
  
  # Remove ranks that are placeholders for missing ranks
  placeholder_ranks <- c('clade', 'no rank')
  table_class_data <- lapply(table_class_data, function(x) x[! x$rank %in% placeholder_ranks, , drop = FALSE])
  
  if (ranks_as_cols) {
    # Reformat with ranks as columns
    parts <- lapply(unique_tax_ids, function(id) {
      values <- table_class_data[[id]][['name']]
      names(values) <- table_class_data[[id]]$rank
      tibble::as_tibble(c(taxon_id = id, as.list(values)))
    })
    
    # Fill in ranks that don't exist in all classifications
    all_cols <- unlist(lapply(parts, colnames))
    col_index <- unlist(lapply(parts, function(x) rev(seq_along(x))))
    unique_cols <- unique(all_cols[order(col_index, decreasing = TRUE)])
    output <- do.call(rbind, lapply(parts, function(part) {
      out <- do.call(data.frame, as.list(rep(NA_character_, length(unique_cols))))
      colnames(out) <- unique_cols
      out[colnames(part)] <- part
      return(out)
    }))
  } else {
    output <- do.call(rbind, table_class_data)
  }
  tibble::as_tibble(output)
}


#' Get taxa predicted by sendsketch results
#'
#' Return which taxa are predicted to be present given sendsketch results
#' associated with directory paths containing `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @inheritParams sendsketch_parsed
#'
#' @return A [tibble::tibble()] with one row per taxon found
#' @family parsers
#'
#' @export
sendsketch_taxa_present <- function(path, ani_thresh = c(species = 95, genus = 90, family = 70),
                                    comp_thresh = c(species = 40, genus = 15, family = 5),
                                    sample_id = NULL, only_best = FALSE, simplify = TRUE) {
  tax_data <- sendsketch_taxonomy_parsed(input = path, ranks_as_cols = FALSE, only_shared = FALSE,
                                         sample_id = sample_id, only_best = only_best, simplify = FALSE)
  sketch_data <- sendsketch_parsed(path = path, sample_id = sample_id, only_best = only_best, simplify = FALSE)
  
  filter_and_extract <- function(rank) {
    subset_ids <- sketch_data$TaxID[sketch_data$ANI > ani_thresh[rank] & sketch_data$Complt > comp_thresh[rank]]
    subset_class_data <- tax_data[tax_data$query_id %in% subset_ids & tax_data$rank == rank, , drop = FALSE]
    subset_class_data[, c('id', 'name', 'rank'), drop = FALSE]
  }
  unique(do.call(rbind, lapply(names(ani_thresh), filter_and_extract)))
}


#' Parse pathogensurveillance output files
#'
#' Parse data for specific pathogensurveillance output files
#'
#' @inheritParams find_ps_paths
#' @param simplify If `TRUE`, attempt to combine all of the output found into a
#'   single object instead of returning a list.
#'
#' @export
find_ps_data <- function(path, target, simplify = FALSE) {

  # Find output paths
  path_data <- find_ps_paths(path = path, target = target, simplify = FALSE, all_metadata = TRUE)
  
  # Check that simplify can be used if specified
  if (simplify) {
    if (length(path_data) != 1) {
      stop(call. = FALSE, 'The `simplify` option is only supported for a single `target`.')
    }
    combiner <- unique(path_data[[1]]$combiner)
    if (is.na(combiner)) {
      stop(call. = FALSE, 'The `simplify` option is not supported for this data type.')
    }
  }
  
  # Get parsers
  unique_parsers <- unique(unlist(lapply(path_data, function(x) x$parser)))
  unique_parsers <- unique_parsers[! is.na(unique_parsers)]
  parsers <- lapply(unique_parsers, function(x) eval(parse(text = x)))
  names(parsers) <- unique_parsers
  
  # Parse data
  output <- lapply(path_data, function(path_data_subset) {
    path_data_subset <- path_data_subset[! is.na(path_data_subset$parser), , drop = FALSE]
    paths <- path_data_subset$path
    out <- lapply(seq_len(nrow(path_data_subset)), function(i) {
      parsers[[path_data_subset$parser[i]]](paths[i])
    })
    names(out) <- paths
    return(out)
  })
  output <- unlist(unname(output), recursive = FALSE)
  
  # Simplify data if possible
  if (simplify) {
    combiner_func <- eval(parse(text = combiner))
    output <- combiner_func(output)
  }
  
  return(output)  
}


#' @keywords internal
combine_tables <- function(tables) {
  # Look up input metadata data to differentiate rows once combined
  path_data <- find_ps_paths(names(tables), target = NULL, simplify = TRUE, must_exist = FALSE)
  
  # Add input metadata columns to each table
  cols_to_add <- colnames(path_data)[! colnames(path_data) %in% c('path', 'target')]
  tables <- lapply(names(tables), function(p) {
    out <- tables[[p]]
    for (col in cols_to_add) {
      out[[col]] <- path_data[[col]][path_data$path == p]
    }
    return(out)
  })
  
  # Combine all of the tables into one
  output <- combine_data_frames(tables)
  
  # Put added columns first
  col_order <- c(cols_to_add, colnames(output)[! colnames(output) %in% cols_to_add])
  output <- output[, col_order, drop = FALSE]
  
  return(output)
}


#' @keywords internal
parse_tsv <- function(path) {
  output <- utils::read.table(path, header = TRUE, sep = '\t', comment.char = '')
  tibble::as_tibble(output)
}


#' @keywords internal
parse_csv <- function(path) {
  output <- utils::read.table(path, header = TRUE, sep = ',', comment.char = '')
  tibble::as_tibble(output)
}


#' @keywords internal
parse_ref_list <- function(path) {
  tibble::tibble(
    reference_id = readLines(path)
  )
}


#' @keywords internal
parse_sendsketch <- function(path) {
  output <- suppressWarnings(utils::read.table(path, sep = '\t', skip = 2, check.names = FALSE, header = TRUE))
  tibble::as_tibble(output)
}


#' @keywords internal
parse_matrix_csv <- function(path) {
  mat <- as.matrix(utils::read.table(path, sep = ',', check.names = FALSE, header = TRUE))
  rownames(mat) <- colnames(mat)
  mat[mat == 0] <- NA
  return(mat)
}

#' @keywords internal
parse_matrix_tsv <- function(path) {
  mat <- as.matrix(utils::read.table(path, sep = '\t', check.names = FALSE, header = TRUE))
  rownames(mat) <- colnames(mat)
  mat[mat == 0] <- NA
  return(mat)
}


#' @keywords internal
parse_tree <- function(path) {
  tree <- ape::read.tree(path)
  tree <- phytools::midpoint_root(tree)
  return(tree)
}


#' @keywords internal
parse_fasta <- function(path) {
  suppressWarnings(ape::read.dna(path, format = "fasta"))
}


#' @keywords internal
parse_run_info <- function(path) {
  run_info <- yaml::read_yaml(path, readLines.warn = FALSE)
  run_info <- lapply(run_info, function(x) {
    if (is.null(x)) {
      return(NA_character_)
    } else {
      return(paste0(x, collapse = ';'))
    }
  })
  tibble::as_tibble(run_info)
}


#' @keywords internal
parse_version_info <- function(path) {
  raw_version_data <- yaml::read_yaml(path)
  tibble::tibble(
    module = rep(names(raw_version_data), vapply(raw_version_data, FUN.VALUE = numeric(1), length)),
    program = unname(unlist(lapply(raw_version_data, names))),
    version = unname(unlist(raw_version_data))
  )
}


#' Find and parse NCBI reference metadata
#'
#' Returns the parsed NCBI reference metadata for all references considered for
#' download.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#'
#' @return A table of reference metadata
#' 
#' @export
parse_ncbi_assembly_metadata <- function(path = NULL) {
  # Read TSV file
  assem_data <- parse_tsv(path)
  
  # Ensure correct column types
  assem_data$tax_id <- as.character(assem_data$tax_id)
  
  # Add taxon info columns
  assem_data$organism_name <- gsub(assem_data$organism_name, pattern = '[', replacement = '', fixed = TRUE)
  assem_data$organism_name <- gsub(assem_data$organism_name, pattern = ']', replacement = '', fixed = TRUE)
  genus_prefixes <- c('Candidatus')
  genus_prefix_pattern <- paste0(genus_prefixes, collapse = '|')
  binomial_pattern <- paste0('^((?:', genus_prefix_pattern, ')?) ?([a-zA-Z0-9.]+) ?([a-zA-Z0-9.]+) ?(.*)$')
  assem_data$species <- trimws(gsub(assem_data$organism_name, pattern = binomial_pattern, replacement = '\\1 \\2 \\3', ignore.case = TRUE))
  assem_data$genus <- trimws(gsub(assem_data$organism_name, pattern = binomial_pattern, replacement = '\\1 \\2', ignore.case = TRUE))
  assem_data$is_refseq <- assem_data$source_database == 'SOURCE_DATABASE_REFSEQ'
  assem_data$is_binomial <- is_latin_binomial(assem_data$species)
  assem_data$assembly_level <- factor(assem_data$assembly_level, ordered = TRUE,
                                      levels = c("Contig", "Scaffold", "Chromosome", "Complete Genome"))
  
  tibble::as_tibble(assem_data)
}



