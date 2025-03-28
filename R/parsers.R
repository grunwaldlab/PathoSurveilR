#' Get parsed pipeline status data
#'
#' Return a [tibble::tibble()] (table) with the status messages produced by the
#' pathogensuriveillance pipeline. The contents of all status message files
#' found in the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with the messages from all input paths
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_parsed(path)
#'
#' @export
status_message_parsed <- function(path) {
  path_data <- status_message_path_data(path)
  if (nrow(path_data) > 0) {
    output <- do.call(rbind, lapply(1:nrow(path_data), function(index) {
      table <- utils::read.table(path_data$path[index], sep = '\t', check.names = FALSE, header = TRUE)
      table$report_group_id <- path_data$report_group_id[index]
      return(table)
    }))
    output <- output[, c("sample_id", "reference_id", "report_group_id", "workflow", "level", "message")]
  } else {
    output <- tibble::tibble(
      sample_id = character(0),
      reference_id = character(0),
      report_group_id = character(0),
      workflow = character(0),
      level = character(0),
      message = character(0)
    )
  }
  return(output)
}


#' Get summary of pipeline status data
#'
#' Return a [tibble::tibble()] (table) with the numbers of issues encountered by the
#' pipeline. The contents of all status message files found in the given paths
#' will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with counts of messages attributes
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_parsed_summary(path)
#'
#' @export
status_message_parsed_summary <- function(path) {
  message_data <- status_message_parsed(path)
  output <- do.call(rbind, lapply(split(message_data, message_data$workflow), function(subset) {
    data.frame(
      workflow = unique(subset$workflow),
      errors = sum(subset$level == 'ERROR'),
      warnings = sum(subset$level == 'WARNING'),
      notes = sum(subset$level == 'NOTE'),
      samples = length(unique(subset$sample_id)),
      references = length(unique(subset$reference_id)),
      report_groups = length(unique(subset$report_group_id))
    )
  }))
  output <- tibble::as_tibble(output)
  return(output)
}


#' Get parsed sample metadata
#'
#' Return a [tibble::tibble()] (table) with sample metadata. The contents of all
#' sample metadata files found in the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with the sample metadata
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_parsed(path)
#'
#' @export
sample_meta_parsed <- function(path) {
  output <- do.call(rbind, lapply(sample_meta_path(path), utils::read.csv, check.names = FALSE, sep = '\t'))
  output[] <- lapply(output, function(col_data) ifelse(col_data == 'null', NA_character_, col_data))
  output <- tibble::as_tibble(output)
  output <- unique(output)
  return(output)
}


#' Get parsed reference metadata
#'
#' Return a [tibble::tibble()] (table) with reference metadata. The contents of
#' all reference metadata files found in the given paths will be combined. This
#' only contains data about references used by at least one step of the pipeline.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [tibble::tibble()] with reference metadata
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_parsed(path)
#'
#' @export
ref_meta_parsed <- function(path) {
  output <- do.call(rbind, lapply(ref_meta_path(path), utils::read.csv, check.names = FALSE, sep = '\t'))
  output[] <- lapply(output, function(col_data) ifelse(col_data == 'null', NA_character_, col_data))
  output <- tibble::as_tibble(output)
  output <- unique(output)
  return(output)
}


#' Get estimated ANI distance matrix
#'
#' Return a list of [base::data.frame()]s with the estimated pairwise ANI values
#' calculated by sourmash.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return a `list` of ANI matrices
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_matrix_parsed(path)
#'
#' @export
estimated_ani_matrix_parsed <- function(path) {
  ani_matrix_paths <- estimated_ani_matrix_path(path)
  output <- lapply(ani_matrix_paths, function(path) {
    ani_matrix <- utils::read.csv(path, check.names = FALSE)
    rownames(ani_matrix) <- colnames(ani_matrix)
    ani_matrix[ani_matrix == 0] <- NA
    return(ani_matrix)
  })
  names(output) <- ani_matrix_paths
  return(output)
}

#' Get POCP matrix
#'
#' Return a list of [base::data.frame()]s with the POCP values based on a core
#' gene analysis using pirate.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return a `list` of POCP matrices
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_matrix_parsed(path)
#'
#' @export
pocp_matrix_parsed <- function(path) {
  matrix_paths <- pocp_matrix_path(path)
  output <- lapply(matrix_paths, function(path) {
    pocp_matrix <- utils::read.csv(path, check.names = FALSE, sep = '\t')
    rownames(pocp_matrix) <- colnames(pocp_matrix)
    return(pocp_matrix)
  })
  names(output) <- matrix_paths
  return(output)
}


#' Get parsed report groups
#'
#' Return a [tibble::tibble()] (table) with report group ID. The groups found in
#' the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return A [base::character()] vector with the groups from all input paths
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' run_info_parsed(path)
#'
#' @export
run_info_parsed <- function(path) {
  path_data <- run_info_path_data(path)
  do.call(rbind, lapply(path_data$path, function(p) {
    run_info <- yaml::read_yaml(p)
    run_info <- lapply(run_info, function(x) {
      if (is.null(x)) {
        return(NA_character_)
      } else {
        return(paste0(x, collapse = ';'))
      }
    })
    tibble::as_tibble(run_info)
  }))
}


#' Get parsed sendsketch results
#'
#' Return a [tibble::tibble()] (table) with the results from sendsketch. The
#' contents of all sendsketch files found in the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param only_best Only return the best hit for each combination of report
#'   group and sample.
#'
#' @return A [tibble::tibble()] with the sendsketch output combined
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_parsed(path)
#'
#' @export
sendsketch_parsed <- function(path, only_best = FALSE) {
  path_data <- sendsketch_path_data(path)

  # If no files are found, return an empty tibble
  if (nrow(path_data) == 0) {
    return(tibble::tibble())
  }

  # Internal function to parse a single file
  parse_one_file <- function(path, report_group_id, sample_id) {
    data <- utils::read.table(path, sep = '\t', skip = 2, check.names = FALSE, header = TRUE)
    numeric_cols <- c(5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                      27, 28, 29, 30, 31, 32, 33)
    data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)
    return(cbind(
      sample_id = rep(sample_id, nrow(data)),
      report_group_id = rep(report_group_id, nrow(data)),
      data
    ))
  }

  # Parse all files and combine them into one data frame
  sketch_data <- do.call(rbind, lapply(seq_len(nrow(path_data)), function(i) {
    parse_one_file(path_data$path[i], path_data$report_group_id[i], path_data$sample_id[i])
  }))

  # Convert percentage fields from character to numeric
  sketch_data$WKID <- as.numeric(gsub("%", "", sketch_data$WKID))
  sketch_data$ANI <- as.numeric(gsub("%", "", sketch_data$ANI))
  sketch_data$Complt <- as.numeric(gsub("%", "", sketch_data$Complt))

  # Filter for top hits
  if (only_best) {
    sketch_data <- sendsketch_best_hits(sketch_data)
  }

  return(tibble::as_tibble(sketch_data))
}


#' Get parsed sendsketch taxonomy
#'
#' Return the taxonomic classificaiton in sendsketch results associated with
#' directory paths containing `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param remove_ranks If `TRUE`, remove the rank information from the taxonomy.
#' @param only_best Only return the best hit for each combination of report
#'   group and sample. 
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
sendsketch_taxonomy_parsed <- function(path, remove_ranks = FALSE, only_best = FALSE) {
  sendsketch_data <- sendsketch_parsed(path, only_best = TRUE)
  classifications <- sendsketch_data$taxonomy
  if (remove_ranks) {
    classifications <- gsub(classifications, pattern = ';?[a-z]+:', replacement = ';')
    classifications <- gsub(classifications, pattern = '^;', replacement = '')
  }
  names(classifications) <- sendsketch_data$sample_id
  classifications <- classifications[! duplicated(paste0(classifications, names(classifications)))]
  return(classifications)
}


#' Get parsed sendsketch taxonomy data
#'
#' Return the taxonomic data in sendsketch results associated with directory
#' paths containing `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param only_best Only return the best hit for each combination of report
#'   group and sample.
#' @param only_shared If `TRUE`, only return the data for that are present in
#'   all of the inputs.
#'
#' @return A [tibble::tibble()] with taxonomy data, with columns corresponding
#'   to ranks.
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_taxonomy_data_parsed(path)
#' sendsketch_taxonomy_data_parsed(path, only_best = TRUE)
#' sendsketch_taxonomy_data_parsed(path, only_shared = TRUE)
#'
#' @export
sendsketch_taxonomy_data_parsed <- function(path, only_best = FALSE, only_shared = FALSE) {
  # Find and parse sendsketch data
  sendsketch_data <- sendsketch_parsed(path, only_best = only_best)

  # Reformat with ranks as columns
  parts <- lapply(1:nrow(sendsketch_data), function(index) {
    split_tax <- strsplit(sendsketch_data$taxonomy[index], split = ';', fixed = TRUE)[[1]]
    taxon_names <- gsub(split_tax, pattern = '^[a-z]+:', replacement = '')
    ranks <- gsub(split_tax, pattern = '^([a-z]*):?.+$', replacement = '\\1')
    names(taxon_names) <- ifelse(ranks == '', 'tip', ranks)
    tibble::as_tibble(as.list(c(sample_id = sendsketch_data$sample_id[index], taxon_names)))
  })

  # Fill in ranks that dont exist in all classifications
  all_cols <- unlist(lapply(parts, colnames))
  col_index <- unlist(lapply(parts, function(x) rev(seq_along(x))))
  unique_cols <- unique(all_cols[order(col_index, decreasing = TRUE)])
  output <- do.call(rbind, lapply(parts, function(part) {
    out <- do.call(data.frame, as.list(rep(NA_character_, length(unique_cols))))
    colnames(out) <- unique_cols
    out[colnames(part)] <- part
    return(out)
  }))

  # Remove columns not shared by all samples
  if (only_shared) {
    col_has_na <- apply(output, MARGIN = 2, function(col) any(is.na(col)))
    output <- output[, ! col_has_na]
  }

  return(tibble::as_tibble(output))
}


#' Parse Software Version Metadata from YAML File
#'
#' Reads a YAML file containing software version information and transforms it
#' into a tibble. Each software module and program, along with its version, is
#' extracted and organized.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return A [tibble::tibble()] with columns `module`, `program`, and `version`
#'   detailing the software versions.
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' software_version_parsed(path)
#'
#' @export
software_version_parsed <- function(path) {
  parse_one <- function(version_path) {
    raw_version_data <- yaml::read_yaml(version_path)
    version_data <- tibble::tibble(
      module = rep(names(raw_version_data), vapply(raw_version_data, FUN.VALUE = numeric(1), length)),
      program = unname(unlist(lapply(raw_version_data, names))),
      version = unname(unlist(raw_version_data))
    )
    return(version_data)
  }
  unique(do.call(rbind, lapply(software_version_path(path), parse_one)))
}


#' Get core gene phylogeny
#'
#' Return a list of [ape::phylo()] objects named by file path by searching for
#' folders containing pathogensurveillance output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_tree_parsed(path)
#'
#' @export
core_tree_parsed <- function(path) {
  tree_parsed(core_tree_path(path))
}


#' Get SNP phylogeny
#'
#' Return a list of [ape::phylo()] objects named by file path by searching for
#' folders containing pathogensurveillance output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param rename If `TRUE`, rename the tip labels to sample IDs and reference IDs.
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_tree_parsed(path)
#'
#' @export
variant_tree_parsed <- function(path, rename = TRUE) {
  path_data <- variant_tree_path_data(path)
  trees <- tree_parsed(path_data$path)

  # Rename tree tips to sample/reference IDs
  if (rename) {
    trees <- lapply(seq_len(length(trees)), function(i) {
      tree <- trees[[i]]
      ref_id <- path_data$ref_id[path_data$path == names(trees)[i]]
      tree$tip.label <- gsub(tree$tip.label, pattern = paste0('^', ref_id, '_'), replacement = '')
      tree$tip.label[tree$tip.label == 'REF'] <- ref_id
      return(tree)
    })
    names(trees) <- path_data$path
  }

  return(trees)
}


#' Get BUSCO phylogeny
#'
#' Return a list of [ape::phylo()] objects named by file path by searching for
#' folders containing pathogensurveillance output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_tree_parsed(path)
#'
#' @export
busco_tree_parsed <- function(path) {
  path_data <- busco_tree_path_data(path)
  trees <- tree_parsed(path_data$path)
  return(trees)
}


#' Get mulitgene phylogenies
#'
#' Return a list of [ape::phylo()] objects named by file path by searching for
#' folders containing pathogensurveillance output. This includes all types of
#' multigene trees produced by pathogensurveillance, such as the core gene
#' phylogeny and the busco gene phylogeny.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' multigene_tree_parsed(path)
#'
#' @export
multigene_tree_parsed <- function(path) {
  return(c(
    busco_tree_parsed(path),
    core_tree_parsed(path)
  ))
}


#' Get phylogenies using a function
#'
#' Return a list of [ape::phylo()] objects named by file path by searching for
#' folders containing pathogensurveillance output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @keywords internal
tree_parsed <- function(path) {
  output <- lapply(path, function(tree_path) {
    tree <- ape::read.tree(tree_path)
    tree <- phytools::midpoint_root(tree)
    return(tree)
  })
  names(output) <- path
  return(output)
}


#' Find and parse SNP alignments
#'
#' Returns parsed SNP alignments for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param rename If `TRUE`, rename the sequence labels as sample IDs and reference IDs.
#' 
#' @return list of [ape::DNAbin()], named by alignment file
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_align_parsed(path)
#'
#' @export
variant_align_parsed <- function(path, rename = TRUE) {
  align_data <- variant_align_path_data(path)
  sample_data <- sample_meta_parsed(path)
  ref_data <- ref_meta_parsed(path)
  output <- lapply(seq_len(nrow(align_data)), function(i) {
    out <- suppressWarnings(ape::read.dna(align_data$path[i], format = "fasta"))
    if (rename && ! is.null(out)) {
      rownames(out) <- gsub(rownames(out), pattern = paste0('^', align_data$ref_id[i], '_'), replacement = '')
      rownames(out)[rownames(out) == 'REF'] <- align_data$ref_id[i]
    }
    return(out)
  })
  names(output) <- align_data$path
  return(output)
}


#' Find and parse NCBI reference metadata
#'
#' Returns the parsed NCBI reference metadata for all references considered for
#' download.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param family The names of families to return data for.
#' @param json_path Path to one or more JSON files containing reference metadata
#'   to parse. It is assumed that these files are named by families. This is an
#'   alternative input to the `path` parameter and both can not be used at once.
#'
#' @return A table of reference metadata
#' 
#' @export
considered_ref_meta_parsed <- function(path = NULL, family = NULL, json_path = NULL) {
  # Check input type parameters
  if (sum(c(is.null(json_path), is.null(path))) != 1) {
    stop('Either the `path` or `json_path` parameters must be used but not both')
  }
  if (! is.null(json_path) && is.null(family)) {
    warning('The `family` parameter has no effect when the `json_path` parameter is used for input.')
  }
  
  # Find data on input files and filter by 'family' argument
  if (is.null(json_path)) {
    path_data <- considered_ref_meta_path_data(path)
    if (is.null(family)) {
      family <- path_data$family
    } else {
      if (! is.character(family)) {
        stop('The argument `family` must be a character vector')
      }
      invalid_families <- family[! family %in% path_data$family]
      if (length(invalid_families) > 0) {
        stop(paste0(
          'The argument `family` contains families not found in the input data:\n  ',
          paste0(invalid_families, collapse = ', '), '\n',
          'Data was found for the following families:\n  ',
          paste0(path_data$family, collapse = ', '), '\n'
        ))
      }
    }
    json_path <- path_data$path[match(family, path_data$family)]
  } else {
    family <- gsub(basename(json_paths), pattern = '.json', replacement = '', fixed = TRUE)
  }
  
  # Parse JSON files into a table
  json_data <- lapply(seq_along(json_path), function(index) {
    parsed_json <- RcppSimdJson::fparse(readLines(json_path[index]), always_list = TRUE)
    output <- do.call(rbind, lapply(parsed_json, function(assem_data) {
      attributes <- assem_data$assembly_info$biosample$attributes
      hosts <- paste0(attributes$value[attributes$name == 'host'], collapse = ';')
      data_parts <- list(
        accession = assem_data$accession,
        assembly_level = assem_data$assembly_info$assembly_level,
        assembly_status = assem_data$assembly_info$assembly_status,
        assembly_type = assem_data$assembly_info$assembly_type,
        hosts = ifelse(hosts == '', NA_character_, hosts),
        organism_name = gsub(assem_data$organism$organism_name, pattern = '\\[|\\]', replacement = ''),
        tax_id = as.character(assem_data$organism$tax_id),
        contig_l50 = as.numeric(assem_data$assembly_stats$contig_l50),
        contig_n50 = as.numeric(assem_data$assembly_stats$contig_n50),
        coverage = as.numeric(sub(assem_data$assembly_stats$genome_coverage, pattern = 'x$', replacement = '')),
        number_of_component_sequences = as.numeric(assem_data$assembly_stats$number_of_component_sequences),
        number_of_contigs = as.numeric(assem_data$assembly_stats$number_of_contigs),
        total_ungapped_length = as.numeric(assem_data$assembly_stats$total_ungapped_length),
        total_sequence_length = as.numeric(assem_data$assembly_stats$total_sequence_length),
        source_database = assem_data$source_database,
        is_type = "type_material" %in% names(assem_data),
        is_annotated = "annotation_info" %in% names(assem_data),
        is_atypical = "atypical" %in% names(assem_data$assembly_info),
        checkm_completeness = assem_data$checkm_info$completeness,
        checkm_contamination = assem_data$checkm_info$contamination
      )
      data_parts[sapply(data_parts, length) == 0 | sapply(data_parts, is.null)] <- NA
      as.data.frame(data_parts)
    }))
    if (!is.null(output)) {
      output$family <- rep(family[index], nrow(output))
      output$genus <- gsub(output$organism_name, pattern = '([a-zA-Z0-9.]+) (.*)', replacement = '\\1')
      output$species <- gsub(output$organism_name, pattern = '([a-zA-Z0-9.]+) ([a-zA-Z0-9.]+) (.*)', replacement = '\\1 \\2')
    }
    return(output)
  })
  assem_data <- do.call(rbind, json_data)
  if (is.null(assem_data)) {
    assem_data <- data.frame(
      accession = character(0),
      assembly_level = character(0),
      assembly_status = character(0),
      assembly_type = character(0),
      hosts = character(0),
      organism_name = character(0),
      tax_id = character(0),
      contig_l50 = numeric(0),
      contig_n50 = numeric(0),
      coverage = numeric(0),
      number_of_component_sequences = numeric(0),
      number_of_contigs = numeric(0),
      total_ungapped_length = numeric(0),
      total_sequence_length = numeric(0),
      source_database = character(0),
      is_type = logical(0),
      is_annotated = logical(0),
      is_atypical = logical(0),
      checkm_completeness = numeric(0),
      checkm_contamination = numeric(0),
      family = character(0),
      genus = character(0),
      species = character(0)
    )
  }
  return(tibble::as_tibble(assem_data))
}
