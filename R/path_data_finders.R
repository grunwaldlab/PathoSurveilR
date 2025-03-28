#' Find sample metadata path data
#'
#' Return the file path data to the TSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_path_data(path)
#'
#' @export
sample_meta_path_data <- function(path) {
  make_path_data_with_group(path, sample_meta_path)
}

#' Find reference metadata path data
#'
#' Return the file path data to the TSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_path_data(path)
#'
#' @export
ref_meta_path_data <- function(path) {
  make_path_data_with_group(path, ref_meta_path)
}

#' Find the BUSCO tree path data
#'
#' Return a table with file paths to the Newick formatted tree produced by
#' comparing BUSCO genes of samples and references for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_tree_path_data(path)
#'
#' @export
busco_tree_path_data <- function(path) {
  make_path_data_with_group(path, busco_tree_path)
}

#' Find the BUSCO analysis reference path data
#'
#' Return a table with the file path to the TSV with a list of references used
#' in the BUSCO analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_ref_path_data(path)
#'
#' @export
busco_ref_path_data <- function(path) {
  make_path_data_with_group(path, busco_ref_path)
}

#' Find the core gene analysis reference path data
#'
#' Return a table with the file path to the TSV with the list of references used
#' in the core gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_ref_path_data(path)
#'
#' @export
core_ref_path_data <- function(path) {
  make_path_data_with_group(path, core_ref_path)
}

#' Find the run info file path data
#'
#' Return a table with the file containing the information about a run in a
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' run_info_path_data(path)
#'
#' @export
run_info_path_data <- function(path) {
  make_path_data_with_group(path, run_info_path)
}

#' Find the variant analysis reference path data
#'
#' Return a table with the file path to the TSV with the list of references used
#' in the variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_ref_path_data(path)
#'
#' @export
variant_ref_path_data <- function(path) {
  make_path_data_with_group(path, variant_ref_path)
}

#' Find the status message TSV path data
#'
#' Return a table with the file path to the TSV with the status reports,
#' warnings, and errors for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_path_data(path)
#'
#' @export
status_message_path_data <- function(path) {
  make_path_data_with_group(path, status_message_path)
}

#' Find the POCP matrix path data
#'
#' Return a table with the file path to the TSV with the POCP (percent of
#' conserved protein matrix for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_matrix_path_data(path)
#'
#' @export
pocp_matrix_path_data <- function(path) {
  make_path_data_with_group(path, pocp_matrix_path)
}

#' Find the estimated ANI matrix path data
#'
#' Return a table with the file path to the CSV with the approximate ANI
#' (average nucleotide identity) matrix estimated by sourmash for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_matrix_path_data(path)
#'
#' @export
estimated_ani_matrix_path_data <- function(path) {
  make_path_data_with_group(path, estimated_ani_matrix_path)
}

#' Find the software version path data
#'
#' Return a table with the file path to the YAML file with the versions of
#' software used for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' software_version_path_data(path)
#'
#' @export
software_version_path_data <- function(path) {
  make_path_data_with_group(path, software_version_path)
}

#' Find the core gene analysis path data
#'
#' Return a table with file paths to the Newick formatted trees produced by the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_tree_path_data(path)
#'
#' @export
core_tree_path_data <- function(path) {
  output <- make_path_data_with_group(path, core_tree_path)
  output$cluster_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(basename(output$path[index]), pattern = paste0('^', output$report_group_id[index], '_cluster_'), replacement = '')
    sub(id, pattern = '\\.treefile$', replacement = '')
  }))
  return(output)
}

#' Find the considered NCBI reference metadata path data
#'
#' Return a table with the file paths of the TSVs of metadata for references
#' considered by the pipeline for download for a given pathogensurveillance
#' output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @export
considered_ref_meta_path_data <- function(path) {
  output <- make_path_data_with_group(path, considered_ref_meta_path)
  output$family <- unlist(lapply(basename(output$path), function(file_name) {
    sub(file_name, pattern = '\\.json$', replacement = '')
  }))
  output <- output[c('report_group_id', 'family', 'path')]
  return(output)
}

#' Find the downloaded reference metadata path data
#'
#' Return a table with the file paths of the TSVs of metadata for references
#' selected and downloaded for each sample for a given pathogensurveillance
#' output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' selected_ref_meta_path_data(path)
#'
#' @export
selected_ref_meta_path_data <- function(path) {
  output <- make_path_data_with_group(path, selected_ref_meta_path)
  output$sample_id <- unlist(lapply(basename(output$path), function(file_name) {
    sub(file_name, pattern = '\\.tsv$', replacement = '')
  }))
  return(output)
}

#' Find the sendsketch result path data
#'
#' Return a table with the file paths of the sendsketch results for each sample
#' for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param sample_id One or more sample IDs to return data for. By default, data
#'   for all samples is returned.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_path_data(path)
#'
#' @export
sendsketch_path_data <- function(path, sample_id = NULL) {
  output <- make_path_data_with_group(path, sendsketch_path)
  output$sample_id <- unlist(lapply(basename(output$path), function(file_name) {
    sub(file_name, pattern = '\\.txt$', replacement = '')
  }))
  output <- output[c('sample_id', 'report_group_id', 'path')]
  
  # Filter by sample ID
  if (! is.null(sample_id)) {
    output <- path_data[output$sample_id %in% sample_id, , drop = FALSE]
  }
  
  return(output)
}

#' Find the SNP alignment path data
#'
#' Return a table with the file paths of the SNP alignments for each reference
#' used in the variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_align_path_data(path)
#'
#' @export
variant_align_path_data <- function(path) {
  output <- make_path_data_with_group(path, variant_align_path)
  file_names <- basename(output$path)
  output$ref_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(file_names[index], pattern = paste0('^', output$report_group_id[index], '_'), replacement = '')
    sub(id, pattern = '\\.fasta$', replacement = '')
  }))
  return(output)
}

#' Find the SNP tree path data
#'
#' Return a table with the file paths of the SNP tree for each reference used in
#' the variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @return `tibble` with `report_group_id` and `path` columns
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_tree_path_data(path)
#'
#' @export
variant_tree_path_data <- function(path) {
  output <- make_path_data_with_group(path, variant_tree_path)
  file_names <- basename(output$path)
  output$ref_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(file_names[index], pattern = paste0('^', output$report_group_id[index], '_'), replacement = '')
    sub(id, pattern = '\\.treefile$', replacement = '')
  }))
  output$cluster_id <- unlist(lapply(1:nrow(output), function(index) {
    id <- sub(file_names[index], pattern = paste0('^', output$report_group_id[index], '_', output$ref_id[index], '_'), replacement = '')
    sub(id, pattern = '\\.treefile$', replacement = '')
  }))
  return(output)
}

#' Make a table with paths and report group
#'
#' @keywords internal
make_path_data_with_group <- function(path, path_func) {
  out_paths <- path_func(path)

  tibble::tibble(
    report_group_id = find_path_report_group(out_paths),
    path = out_paths
  )
}

#' Determine the report group for files
#'
#' @keywords internal
find_path_report_group <- function(path) {
  find_one <- function(p) {
    while (p != dirname(p)) {
      if ("pathogensurveillance_run_info.yml" %in% list.files(p)) {
        return(yaml::read_yaml(file.path(p, "pathogensurveillance_run_info.yml"))$group_id)
      }
      p <- dirname(p)
    }
    return(NA_character_)
  }
  if (length(path) == 0) {
    return(character(0))
  } else {
    return(unlist(lapply(path, find_one)))
  }
}
