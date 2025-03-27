#' Find sample metadata file path
#'
#' Return the file path to the TSV containing the sample metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param from_report_input If `TRUE` (the default), return the version of the
#'   metadata passed to the report rather than the first parsed version of the
#'   user metadata used for all report groups.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_path(path)
#'
#' @export
sample_meta_path <- function(path, from_report_input = TRUE) {
  if (from_report_input) {
    return(find_static_path(path, output_path = character(0), report_path = 'sample_data.tsv'))
  } else {
    return(find_static_path(path, output_path = 'pipeline_info/sample_metadata.tsv', report_path = character(0)))
  }
}

#' Find reference metadata file path
#'
#' Return the file path to the TSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @param from_report_input If `TRUE` (the default), return the version of the
#'   metadata passed to the report rather than the first parsed version of the
#'   user metadata used for all report groups.
#' @return character vector of length 1
#' @family path finders
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_path(path)
#'
#' @export
ref_meta_path <- function(path, from_report_input = TRUE) {
  if (from_report_input) {
    return(find_static_path(path, output_path = character(0), report_path = 'reference_data.tsv'))
  } else {
    return(find_static_path(path, output_path = 'pipeline_info/reference_metadata.tsv', report_path = character(0)))
  }
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
  find_static_path(path, output_path = 'iqtree2_busco', report_path = 'busco_trees')
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
  find_static_path(path, output_path = 'assign_busco_references', report_path = 'busco_tree_references.tsv')
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
  find_static_path(path, output_path = 'assign_core_references', report_path = 'core_gene_tree_references.tsv')
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
  find_static_path(path, output_path = NULL, report_path = 'pathogensurveillance_run_info.yml')
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
  find_static_path(path, output_path = 'assign_mapping_reference', report_path = 'mapping_references.tsv')
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
  find_static_path(path, output_path = NULL, report_path = 'messages.tsv')
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
  find_static_path(path, output_path = 'calculate_pocp', report_path = 'pocp.tsv')
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
  output <- find_static_path(path, output_path = 'sourmash_compare', report_path = 'sourmash_ani_matrix.csv')
  output[endsWith(output, '.csv')]
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
  find_static_path(path, output_path = NULL, report_path = 'versions.yml')
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
  find_static_path(path, output_path = 'iqtree2_core', report_path = 'core_gene_trees')
}

#' Find the considered NCBI reference metadata paths
#'
#' Return the file paths of the JSONs of metadata for references considered by
#' the pipeline for download for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @return character vector of length 1
#' @family path finders
#'
#' @export
considered_ref_meta_path <- function(path) {
  find_static_path(path, output_path = 'find_assemblies', report_path = 'ncbi_reference_data')
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
  find_static_path(path, output_path = 'pick_assemblies', report_path = 'selected_references')
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
  find_static_path(path, output_path = 'bbmap_sendsketch', report_path = 'sendsketch')
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
  find_static_path(path, output_path = 'vcf_to_snp_align', report_path = 'snp_alignments')
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
  find_static_path(path, output_path = 'iqtree2_snp', report_path = 'snp_trees')
}


#' Find files in pathogensurveillance output
#'
#' Find files in the output of pathogensurveillance, either in the entire output
#' folder or the report input folder.
#'
#' @param path One or more paths to search
#' @param output_path The path of the file/directory relative to the entire
#'   pathogensurveillance output.
#' @param report_path The path of the file/directory relative to the entire
#'   report input directory.
#' @param return_contents If `TRUE` (the default), if `path` is a directory, return the files
#'   in that directory rather than the directory itself.
#'
#' @keywords internal
find_static_path <- function(path, output_path, report_path, return_contents = TRUE, ...) {

  # Look for directories in both the output folder and report input directories
  report_dir_paths <- find_report_input_result_path(path)
  output_dir_paths <- find_output_directory_path(path)
  
  # List files in target directories
  find_one <- function(p, name) {
    expected_path <- file.path(p, name)
    if (file.exists(expected_path)) {
      if (return_contents && dir.exists(expected_path)) {
        out_paths <- list.files(expected_path, full.names = TRUE, ...)
        out_paths <- out_paths[file.exists(out_paths)]
      } else {
        out_paths <- expected_path
      }
    } else {
      out_paths <- character(0)
    }
    return(out_paths)
  }
  if (length(output_path) == 0) {
    output_files <- character(0)
  } else {
    output_files <- lapply(output_dir_paths, function(x) find_one(x, output_path))
  } 
  if (length(report_path) == 0) {
    report_files <- character(0)
  } else {
    report_files <- lapply(report_dir_paths, function(x) find_one(x, report_path))
    if (length(output_files) > 0) { # Ignore files in report input if also present in entire output
      is_redundant <- vapply(seq_len(length(report_dir_paths)), FUN.VALUE = logical(1), function(i) {
        is_nested <- startsWith(report_dir_paths[i], prefix = output_dir_paths)
        if (all(! is_nested)) {
          return(FALSE)
        } 
        return(length(report_files[[i]]) > 0 & length(output_files[[which(is_nested)]]) > 0)
      })
      report_files <- report_files[! is_redundant]
    }
  }
  
  # Return all files found in all directories
  output <- unlist(c(output_files, report_files))
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

