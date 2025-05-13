#' Find sample metadata file path
#'
#' Return the file path to the TSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_path(path)
#'
#' @export
sample_meta_path <- function(path, simplify = FALSE) {
  user_version <- find_path(
    path,
    out_dir_subpath = 'metadata',
    out_dir_pattern = '^sample_metadata\\.tsv$',
    report_dir_subpath = NULL,
    report_dir_pattern =  NULL,
    simplify = FALSE
  )
  report_version <- find_path(
    path,
    out_dir_subpath = 'report_group_data',
    out_dir_pattern = '^.+_inputs/sample_data\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern =  '^sample_data\\.tsv$',
    simplify = FALSE
  )
  has_report_version <- vapply(report_version, FUN.VALUE = numeric(1), length) > 0
  has_user_version <- vapply(user_version, FUN.VALUE = numeric(1), length) > 0
  output <- c(report_version[has_report_version], user_version[has_user_version & ! has_report_version])
  
  if (simplify) {
    output <- unname(unlist(output))
  }
  
  return(output)
}


#' Find reference metadata file path
#'
#' Return the file path to the TSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_path(path)
#'
#' @export
ref_meta_path <- function(path, simplify = FALSE) {
  user_version <- find_path(
    path,
    out_dir_subpath = 'metadata',
    out_dir_pattern = '^reference_metadata\\.tsv$',
    report_dir_subpath = NULL,
    report_dir_pattern =  NULL,
    simplify = simplify
  )
  report_version <- find_path(
    path,
    out_dir_subpath = 'report_group_data',
    out_dir_pattern = '^.+_inputs/reference_data\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern =  '^reference_data\\.tsv$',
    simplify = simplify
  )
  has_report_version <- vapply(report_version, FUN.VALUE = numeric(1), length) > 0
  has_user_version <- vapply(user_version, FUN.VALUE = numeric(1), length) > 0
  output <- c(report_version[has_report_version], user_version[has_user_version & ! has_report_version])
  
  if (simplify) {
    output <- unname(unlist(output))
  }
  
  return(output)
}


#' Find the BUSCO tree file path
#'
#' Return the file path to the Newick formatted tree produced by comparing BUSCO
#' genes of samples and references for a given pathogensurveillance output
#' folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_tree_path(path)
#'
#' @export
busco_tree_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'trees/busco',
    out_dir_pattern = '\\.treefile$',
    report_dir_subpath = 'busco_trees',
    report_dir_pattern = '\\.treefile$',
    simplify = simplify
  )
}


#' Find the BUSCO analysis reference data file path
#'
#' Return the file path to the TSV with a list of references used in the BUSCO
#' analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_ref_path(path)
#'
#' @export
busco_ref_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'reference_data/selected',
    out_dir_pattern = '_busco_references\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern = '^busco_tree_references\\.tsv$',
    simplify = simplify
  )
}


#' Find the core gene analysis reference data file path
#'
#' Return the file path to the TSV with the list of references used in the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_ref_path(path)
#'
#' @export
core_ref_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'reference_data/selected',
    out_dir_pattern = '_core_references\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern = '^core_gene_tree_references\\.tsv$',
    simplify = simplify
  )
}


#' Find the report group file path
#'
#' Return the file containing the name of the report group for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' run_info_path(path)
#'
#' @export
run_info_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'pipeline_info',
    out_dir_pattern = '^pathogensurveillance_run_info\\.yml$',
    report_dir_subpath = '',
    report_dir_pattern = '^pathogensurveillance_run_info\\.yml$',
    simplify = simplify
  )
}


#' Find the variant analysis reference data file path
#'
#' Return the file path to the TSV with the list of references used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_ref_path(path)
#'
#' @export
variant_ref_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'reference_data/selected',
    out_dir_pattern = '_mapping_references\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern = '^mapping_references\\.tsv$',
    simplify = simplify
  )
}


#' Find the status message TSV file path
#'
#' Return the file path to the TSV with the status reports, warnings, and errors
#' for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_path(path)
#'
#' @export
status_message_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'pipeline_info',
    out_dir_pattern = '^messages\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern = '^messages\\.tsv$',
    simplify = simplify
  )
}


#' Find the POCP matrix file path
#'
#' Return the file path to the TSV with the POCP (percent of conserved protein
#' matrix for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_matrix_path(path)
#'
#' @export
pocp_matrix_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'pocp',
    out_dir_pattern = '_pocp\\.tsv$',
    report_dir_subpath = '',
    report_dir_pattern = '^pocp\\.tsv$',
    simplify = simplify
  )
}


#' Find the estimated ANI matrix file path
#'
#' Return the file path to the CSV with the approximate ANI (average nucleotide
#' identity) matrix estimated by sourmash for a given pathogensurveillance
#' output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @param format One of `csv`, `npy`, or `npy.labels.txt`. Formats other than
#'   `csv` are only available when using whole output directories as input.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_matrix_path(path)
#'
#' @export
estimated_ani_matrix_path <- function(path, simplify = FALSE, format = 'csv') {
  find_path(
    path,
    out_dir_subpath = 'sketch_comparisons/ani_matricies',
    out_dir_pattern = paste0('\\.', format, '$'),
    report_dir_subpath = '',
    report_dir_pattern =  paste0('sourmash_ani_matrix\\.', format, '$'),
    simplify = simplify
  )
}


#' Find the software version file path
#'
#' Return the file path to the YAML file with the versions of software used for
#' a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' software_version_path(path)
#'
#' @export
software_version_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'report_group_data',
    out_dir_pattern = '^.+_inputs/versions.yml$',
    report_dir_subpath = '',
    report_dir_pattern = '^versions.yml',
    simplify = simplify
  )
}


#' Find the core gene analysis tree paths
#'
#' Return the file paths to the Newick formatted trees produced by the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_tree_path(path)
#'
#' @export
core_tree_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'trees/core',
    out_dir_pattern = '\\.treefile$',
    report_dir_subpath = 'core_gene_trees',
    report_dir_pattern = '\\.treefile$',
    simplify = simplify
  )
}


#' Find the considered NCBI reference metadata paths
#'
#' Return the file paths of the JSONs of metadata for references considered by
#' the pipeline for download for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @param format One of `json` or `tsv`.
#' @return character vector of length 1
#' @family path finders
#'
#' @export
considered_ref_meta_path <- function(path, format = 'tsv', simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'reference_data/considered',
    out_dir_pattern = paste0('\\.', format, '$'),
    report_dir_subpath = 'ncbi_reference_data',
    report_dir_pattern = paste0('\\.', format, '$'),
    simplify = simplify
  )
}


#' Find the downloaded reference metadata paths
#'
#' Return the file paths of the TSVs of metadata for references selected and
#' downloaded for each sample for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' downloaded_ref_meta_path(path)
#'
#' @export
downloaded_ref_meta_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'reference_data/downloaded',
    out_dir_pattern = '\\.tsv$',
    report_dir_subpath = 'selected_references',
    report_dir_pattern = '\\.tsv$',
    simplify = simplify
  )
}


#' Find the sendsketch result paths
#'
#' Return the file paths of the sendsketch results for each sample for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_path(path)
#'
#' @export
sendsketch_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'sendsketch',
    out_dir_pattern = '\\.txt$',
    report_dir_subpath = 'sendsketch',
    report_dir_pattern = '\\.txt$',
    simplify = simplify
  )
}


#' Find the SNP alignment paths
#'
#' Return the file paths of the SNP alignments for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_align_path(path)
#'
#' @export
variant_align_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'variants',
    out_dir_pattern = '\\.fasta$',
    report_dir_subpath = 'snp_alignments',
    report_dir_pattern = '\\.fasta$',
    simplify = simplify
  )
}


#' Find the SNP tree paths
#'
#' Return the file paths of the SNP tree for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_tree_path(path)
#'
#' @export
variant_tree_path <- function(path, simplify = FALSE) {
  find_path(
    path,
    out_dir_subpath = 'trees/snp',
    out_dir_pattern = '\\.treefile$',
    report_dir_subpath = 'snp_trees',
    report_dir_pattern = '\\.treefile$',
    simplify = simplify
  )
}


#' Search for a specific type of file in pathogensurveillance output
#'
#' Search for a specific type of file in pathogensurveillance output
#'
#' @param path One or more paths to search for pathogensurveillance output,
#'   either entire output directories or report input directories.
#' @param out_dir_subpath A path relative to `path` that does not change. Used
#'   to lessen the number of files that need to be searched through.
#' @param out_dir_pattern A regex matching the path of the file/directory
#'   relative to `path`/`out_dir_subpath`.
#' @param report_dir_subpath A path relative to `path` that does not change.
#'   Used to lessen the number of files that need to be searched through.
#' @param report_dir_pattern A regex matching the path of the file/directory
#'   relative to `path`/`report_dir_subpath`.
#' @param simplify If `TRUE`, return a single vector with all paths rather than
#'   a list of vectors for each output directory.
#'
#' @keywords internal
find_path <- function(path, out_dir_subpath, out_dir_pattern, report_dir_subpath, report_dir_pattern, simplify = TRUE) {
  # Remove trailing slash if present
  path <- sub(path, pattern = '/$', replacement = '')
  
  # Find outputs of the pathogensurveillance pipeline
  output_dir_paths <- find_output_directory_path(path)
  report_dir_paths <- find_report_input_result_path(path)
  
  # Remove any report input paths that are present inside entire output directories
  if (length(output_dir_paths) > 0) {
    report_dir_in_output_dir <- vapply(report_dir_paths, FUN.VALUE = logical(1), function(p) {
      any(startsWith(p, paste0(output_dir_paths, '/')))
    })
    report_dir_paths <- report_dir_paths[! report_dir_in_output_dir]
  }
  
  # Find files of interest inside pathogensurveillance output directories
  output <- list()
  output_names <- character(0)
  if (! (is.null(out_dir_subpath) && is.null(out_dir_pattern)) ) {
    output <- c(output, lapply(output_dir_paths, find_path_in_output_directory,
                               subpath = out_dir_subpath, pattern = out_dir_pattern))
    output_names <- c(output_names, output_dir_paths)
  }
  if (! (is.null(report_dir_subpath) && is.null(report_dir_pattern)) ) {
    output <- c(output, lapply(report_dir_paths, find_path_in_output_directory,
                               subpath = report_dir_subpath, pattern = report_dir_pattern))
    output_names <- c(output_names, report_dir_paths)
  }
  
  # Collapse results to a single vector if specified
  if (simplify) {
    output <- unlist(output)
  } else {
    names(output) <- output_names
  }
  
  return(output)
}


#' Find files in pathogensurveillance output
#'
#' Find files in the output of pathogensurveillance.
#'
#' @param path One or more paths to search.
#' @param subpath A path relative to `path` that does not change. Used to lessen
#'   the number of files that need to be searched through.
#' @param pattern A regex matching the path of the file/directory relative to
#'   `path`/`subpath`.
#'
#' @keywords internal
find_path_in_output_directory <- function(path, subpath, pattern) {
  if (is.null(subpath) || nchar(subpath) == 0) {
    search_path <- path
  } else {
    search_path <- file.path(path, subpath)
  }
  unlist(lapply(search_path, function(p) {
    subsubpaths <- list.files(p, recursive = TRUE, include.dirs = TRUE, all.files = TRUE)
    subsubpaths <- subsubpaths[grepl(subsubpaths, pattern = pattern)]
    out <- file.path(p, subsubpaths)
    out[file.exists(out)]
  }))
}


#' Find subfolders with report group results
#'
#' For one or more paths, return the paths of subfolders (or the given folders)
#' that are the report group output of a pathogensurviellance run.
#'
#' @keywords internal
find_report_input_result_path <- function(path) {
  subdir_paths <- list.dirs(path)
  is_valid_dir <- function(p) {
    group_id_path <- file.path(p, 'pathogensurveillance_run_info.yml')
    return(file.exists(group_id_path))
  }
  subdir_paths <- subdir_paths[unlist(lapply(subdir_paths, is_valid_dir))]
  return(unique(subdir_paths))
}


#' Find pipeline output directory paths
#'
#' For one or more paths, return the paths of subfolders (or the given folders)
#' that are the entire output of a pathogensurviellance run.
#'
#' @keywords internal
find_output_directory_path <- function(path) {
  subdir_paths <- list.dirs(path)
  is_valid_dir <- function(p) {
    pipeline_yml_path <- file.path(p, 'pipeline_info', 'nf_core_pipeline_software_mqc_versions.yml')
    file.exists(pipeline_yml_path) && any(grepl(readLines(pipeline_yml_path), pattern = 'pathogensurveillance'))
  }
  subdir_paths <- subdir_paths[unlist(lapply(subdir_paths, is_valid_dir))]
  return(unique(subdir_paths))
}

