#' Find sample metadata file path
#'
#' Return the file path to the TSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_path(path)
#'
#' @export
sample_meta_path <- function(path) {
  find_static_file_path(path, 'sample_data.tsv', 'sample metadata file')
}

#' Find reference metadata file path
#'
#' Return the file path to the TSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_path(path)
#'
#' @export
ref_meta_path <- function(path) {
  find_static_file_path(path, 'reference_data.tsv', 'reference metadata file')
}

#' Find the BUSCO tree file path
#'
#' Return the file path to the Newick formatted tree produced by comparing BUSCO
#' genes of samples and references for a given pathogensurveillance output
#' folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_tree_path(path)
#'
#' @export
busco_tree_path <- function(path) {
  find_static_dir_path(path, 'busco_trees', 'busco gene tree', file_required = FALSE, dir_required = FALSE)
}

#' Find the BUSCO analysis reference data file path
#'
#' Return the file path to the TSV with a list of references used in the BUSCO
#' analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_ref_path(path)
#'
#' @export
busco_ref_path <- function(path) {
  find_static_file_path(path, 'busco_tree_references.tsv', 'busco anaylsis reference file', file_required = FALSE)
}

#' Find the core gene analysis reference data file path
#'
#' Return the file path to the TSV with the list of references used in the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_ref_path(path)
#'
#' @export
core_ref_path <- function(path) {
  find_static_file_path(path, 'core_gene_tree_references.tsv', 'core gene analysis reference file', file_required = FALSE)
}

#' Find the report group file path
#'
#' Return the file containing the name of the report group for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' run_info_path(path)
#'
#' @export
run_info_path <- function(path) {
  find_static_file_path(path, 'pathogensurveillance_run_info.yml', 'report group name file')
}

#' Find the variant analysis reference data file path
#'
#' Return the file path to the TSV with the list of references used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_ref_path(path)
#'
#' @export
variant_ref_path <- function(path) {
  find_static_file_path(path, 'mapping_references.tsv', 'variant analysis reference file')
}

#' Find the status message TSV file path
#'
#' Return the file path to the TSV with the status reports, warnings, and errors
#' for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_path(path)
#'
#' @export
status_message_path <- function(path) {
  find_static_file_path(path, 'messages.tsv', 'status message data file', file_required = FALSE)
}

#' Find the POCP matrix file path
#'
#' Return the file path to the TSV with the POCP (percent of conserved protein
#' matrix for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_matrix_path(path)
#'
#' @export
pocp_matrix_path <- function(path) {
  find_static_file_path(path, 'pocp.tsv', 'pocp matrix file', file_required = FALSE)
}

#' Find the estimated ANI matrix file path
#'
#' Return the file path to the CSV with the approximate ANI (average nucleotide identity)
#' matrix estimated by sourmash for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_matrix_path(path)
#'
#' @export
estimated_ani_matrix_path <- function(path) {
  find_static_file_path(path, 'sourmash_ani_matrix.csv', 'estimated ANI matrix file')
}

#' Find the software version file path
#'
#' Return the file path to the YAML file with the versions of software used for
#' a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' software_version_path(path)
#'
#' @export
software_version_path <- function(path) {
  find_static_file_path(path, 'versions.yml', 'software version file')
}

#' Find the core gene analysis tree paths
#'
#' Return the file paths to the Newick formatted trees produced by the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_tree_path(path)
#'
#' @export
core_tree_path <- function(path) {
  find_static_dir_path(path, 'core_gene_trees', 'core gene tree', file_required = FALSE, dir_required = FALSE)
}

#' Find the considered NCBI reference metadata paths
#'
#' Return the file paths of the TSVs of metadata for references considered by
#' the pipeline for download for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @export
considered_ref_meta_path <- function(path) {
  find_static_dir_path(path, 'ncbi_reference_data', 'considered NCBI reference metadata')
}

#' Find the downloaded reference metadata paths
#'
#' Return the file paths of the TSVs of metadata for references selected and
#' downloaded for each sample for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' selected_ref_meta_path(path)
#'
#' @export
selected_ref_meta_path <- function(path) {
  find_static_dir_path(path, 'selected_references', 'selected NCBI reference metadata')
}

#' Find the sendsketch result paths
#'
#' Return the file paths of the sendsketch results for each sample for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_path(path)
#'
#' @export
sendsketch_path <- function(path) {
  find_static_dir_path(path, 'sendsketch', 'sendsketch result')
}

#' Find the SNP alignment paths
#'
#' Return the file paths of the SNP alignments for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_align_path(path)
#'
#' @export
variant_align_path <- function(path) {
  find_static_dir_path(path, 'snp_alignments', 'SNP alignment', dir_required = FALSE, file_required = FALSE)
}

#' Find the SNP tree paths
#'
#' Return the file paths of the SNP tree for each reference used in the
#' variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_tree_path(path)
#'
#' @export
variant_tree_path <- function(path) {
  find_static_dir_path(path, 'snp_trees', 'SNP tree', dir_required = FALSE, file_required = FALSE)
}

#' @keywords internal
find_static_file_path <- function(path, file_name, file_description, file_required = TRUE) {

  group_result_paths <- find_report_input_result_path(path)

  find_one <- function(p) {
    output_path <- file.path(p, file_name)
    if (! file.exists(output_path)) {
      if (file_required) {
        stop(
          call. = FALSE,
          'Cannot locate ', file_description, '. It should be located at "',
          output_path, '". Verify that "', path,
          '" is an pathogensurveillance output folder.'
        )
      } else {
        return(character(0))
      }
    }
    return(output_path)
  }
  output <- unlist(lapply(group_result_paths, find_one))
  if (is.null(output)) {
    output <- character(0)
  }
  return(output)
}

#' @keywords internal
find_static_dir_path <- function(path, dir_name, file_description, dir_required = TRUE, file_required = TRUE, ...) {

  group_result_paths <- find_report_input_result_path(path)

  # List files in target directories
  find_one <- function(p) {
    # Check that directory exists
    dir_path <- file.path(p, dir_name)
    if (! file.exists(dir_path)) {
      if (dir_required) {
        stop(
          call. = FALSE,
          'Cannot locate the ', file_description, ' directory. It should be located at "',
          dir_path, '". Verify that "', path,
          '" is an pathogensurveillance output folder.'
        )
      } else {
        return(character(0))
      }
    }
    # List files inside directory
    out_paths <- list.files(dir_path, full.names = TRUE, ...)
    if (file_required && length(out_paths) == 0) {
      stop(
        call. = FALSE,
        'No files found in the ', file_description, ' directory. It should be located at "',
        dir_path, '". Verify that "', path,
        '" is an pathogensurveillance output folder.'
      )
    }
    return(out_paths)
  }
  output <- unlist(lapply(group_result_paths, find_one))
  if (is.null(output)) {
    output <- character(0)
  }
  return(output)
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
  return(subdir_paths)
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
  return(subdir_paths)
}

