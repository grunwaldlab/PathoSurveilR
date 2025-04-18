#' Find sample metadata path data
#'
#' Return the file path data to the TSV containing the sample metadata for a
#' given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_path_data(path)
#'
#' @export
sample_meta_path_data <- function(path, simplify = TRUE) {
  samp_meta_paths <- sample_meta_path(path)
  output <- lapply(samp_meta_paths, function(p) {
    tibble::tibble(
      report_id = sub(p, pattern = '^.*/(.+)_inputs/sample_data\\.tsv$', replacement = '\\1'),
      path = p
    )
  })
  postprocess_path_data(output, simplify = simplify)
}


#' Find reference metadata path data
#'
#' Return the file path data to the TSV containing the reference metadata for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_path_data(path)
#'
#' @export
ref_meta_path_data <- function(path, simplify = TRUE) {
  ref_meta_paths <- ref_meta_path(path)
  output <- lapply(ref_meta_paths, function(p) {
    tibble::tibble(
      report_id = sub(p, pattern = '^.*/(.+)_inputs/reference_data\\.tsv$', replacement = '\\1'),
      path = p
    )
  })
  postprocess_path_data(output, simplify = simplify)
}


#' Find the BUSCO tree path data
#'
#' Return a table with file paths to the Newick formatted tree produced by
#' comparing BUSCO genes of samples and references for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_tree_path_data(path)
#'
#' @export
busco_tree_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = busco_tree_path,
    regex = "^(.+?)_cluster_[0-9]+\\.treefile$",
    id_types = c("report_id")
  )
  output <- lapply(output, function(x) {
    x$cluster_id <- vapply(seq_len(nrow(x)), FUN.VALUE = character(1), function(i) {
      sub(basename(x$path[i]), pattern = paste0('^', x$report_id[i], '_(cluster_[0-9]+)\\.treefile'), replacement = '\\1')
    })
    x[, c('report_id', 'cluster_id', 'path')]
  })
  postprocess_path_data(output, simplify = simplify)
}


#' Find the BUSCO analysis reference path data
#'
#' Return a table with the file path to the TSV with a list of references used
#' in the BUSCO analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_ref_path_data(path)
#'
#' @export
busco_ref_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = busco_ref_path,
    regex = "^(.+?)_busco_references\\.tsv$",
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the core gene analysis reference path data
#'
#' Return a table with the file path to the TSV with the list of references used
#' in the core gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_ref_path_data(path)
#'
#' @export
core_ref_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = core_ref_path,
    regex = "^(.+?)_core_references\\.tsv$",
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the run info file path data
#'
#' Return a table with the file containing the information about a run in a
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' run_info_path_data(path)
#'
#' @export
run_info_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = run_info_path,
    regex = NULL,
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the variant analysis reference path data
#'
#' Return a table with the file path to the TSV with the list of references used
#' in the variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_ref_path_data(path)
#'
#' @export
variant_ref_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = variant_ref_path,
    regex = "^(.+?)_mapping_references\\.tsv$",
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the status message TSV path data
#'
#' Return a table with the file path to the TSV with the status reports,
#' warnings, and errors for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_path_data(path)
#'
#' @export
status_message_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = status_message_path,
    regex = "^(.+?)\\.tsv$",
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the POCP matrix path data
#'
#' Return a table with the file path to the TSV with the POCP (percent of
#' conserved protein matrix for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_matrix_path_data(path)
#'
#' @export
pocp_matrix_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = pocp_matrix_path,
    regex = "^(.+?)_pocp\\.tsv$",
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the estimated ANI matrix path data
#'
#' Return a table with the file path to the CSV with the approximate ANI
#' (average nucleotide identity) matrix estimated by sourmash for a given
#' pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_matrix_path_data(path)
#'
#' @export
estimated_ani_matrix_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = estimated_ani_matrix_path,
    regex = "^(.+?)_comp\\.csv$",
    id_types = c("report_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the software version path data
#'
#' Return a table with the file path to the YAML file with the versions of
#' software used for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' software_version_path_data(path)
#'
#' @export
software_version_path_data <- function(path, simplify = TRUE) {
  version_meta_paths <- software_version_path(path)
  output <- lapply(version_meta_paths, function(p) {
    tibble::tibble(
      report_id = sub(p, pattern = '^.*/(.+)_inputs/versions\\.yml$', replacement = '\\1'),
      path = p
    )
  })
  postprocess_path_data(output, simplify = simplify)
}


#' Find the core gene analysis path data
#'
#' Return a table with file paths to the Newick formatted trees produced by the core
#' gene analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_tree_path_data(path)
#'
#' @export
core_tree_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = core_tree_path,
    regex = "^(.+?)_cluster_[0-9]+\\.treefile$",
    id_types = c("report_id")
  )
  output <- lapply(output, function(x) {
    x$cluster_id <- vapply(seq_len(nrow(x)), FUN.VALUE = character(1), function(i) {
      sub(basename(x$path[i]), pattern = paste0('^', x$report_id[i], '_(cluster_[0-9]+)\\.treefile'), replacement = '\\1')
    })
    x[, c('report_id', 'cluster_id', 'path')]
  })
  postprocess_path_data(output, simplify = simplify)
}


#' Find the considered NCBI reference metadata path data
#'
#' Return a table with the file paths of the TSVs of metadata for references
#' considered by the pipeline for download for a given pathogensurveillance
#' output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @export
considered_ref_meta_path_data <- function(path, simplify = TRUE) {
  ref_meta_paths <- considered_ref_meta_path(path)
  output <- lapply(ref_meta_paths, function(p) {
    tibble::tibble(
      taxon_id = sub(basename(p), pattern = '^(.+?)\\.tsv', replacement = '\\1'),
      path = p
    )
  })
  postprocess_path_data(output, simplify = simplify)
}


#' Find the downloaded reference metadata path data
#'
#' Return a table with the file paths of the TSVs of metadata for references
#' selected and downloaded for each sample for a given pathogensurveillance
#' output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' selected_ref_meta_path_data(path)
#'
#' @export
downloaded_ref_meta_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = downloaded_ref_meta_path,
    regex = "^(.+?)\\.tsv$",
    id_types = c("sample_id")
  )
  postprocess_path_data(output, simplify = simplify)
}

#' Find the sendsketch result path data
#'
#' Return a table with the file paths of the sendsketch results for each sample
#' for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_path_data(path)
#'
#' @export
sendsketch_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = sendsketch_path,
    regex = "^(.+?)\\.txt$",
    id_types = c("sample_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the SNP alignment path data
#'
#' Return a table with the file paths of the SNP alignments for each reference
#' used in the variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_align_path_data(path)
#'
#' @export
variant_align_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = variant_align_path,
    regex = "^(.+?)\\.fasta$",
    id_types = c("report_id", "reference_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Find the SNP tree path data
#'
#' Return a table with the file paths of the SNP tree for each reference used in
#' the variant analysis for a given pathogensurveillance output folder.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_path_data
#' @return A [tibble::tibble()] with one row per path if `simplify = TRUE` or a list of
#'   such [tibble::tibble()]s for each output directory found if `simplify = FALSE`.
#' @family path tables
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_tree_path_data(path)
#'
#' @export
variant_tree_path_data <- function(path, simplify = TRUE) {
  output <- make_path_data_with_group(
    path,
    path_func = variant_tree_path,
    regex = "^(.+?)\\.treefile$",
    id_types = c("report_id", "reference_id")
  )
  postprocess_path_data(output, simplify = simplify)
}


#' Prepare list of tables for output
#'
#' Optionally combine a list of tables into a single table for use with path
#' data finder functions. Also can add columns for which output folders paths
#' came from.
#'
#' @param table_list The list to combine
#' @param simplify If `FALSE` a list of [tibble::tibble()]s are returned named by the
#'   output folder the data was found in. If `TRUE`, all data is combined into a
#'   single [tibble::tibble()].
#'
#' @keywords internal
postprocess_path_data <- function(table_list, simplify) {
  for (outdir_path in names(table_list)) {
    table_list[[outdir_path]]$outdir_path <- outdir_path
  }
  if (simplify) {
    table_list <- do.call(combine_data_frames, table_list)
  }
  return(table_list)
}


#' Make a table with paths and report group
#'
#' @inheritParams find_path
#' @param path_func A path finder function that returns a list of paths
#' @param regex A regular expression with a single capture group that matches
#'   the portion of the file name that is non-constant.
#' @param id_types One or more of `'sample_id'`, `'reference_id'`,
#'   `'report_id'` in the order that they appear in file name.
#' @param sep The separator that is used between elements in the file name.
#'
#' @keywords internal
make_path_data_with_group <- function(path, path_func, regex, id_types, sep = '_') {
  path_data <- path_func(path)
  sample_meta <- sample_meta_parsed(path, simplify = FALSE)
  reference_meta <- ref_meta_parsed(path, simplify = FALSE)
  
  output <- lapply(names(path_data), function(outdir) {
    if (length(path_data[[outdir]]) == 0) {
      # Return empty table if no paths
      out <- data.frame(matrix(vector(), 0, length(id_types) + 1,
                               dimnames = list(c(), c(id_types, 'path'))),
                        stringsAsFactors = FALSE)
    } else {
      if (is.null(regex)) {
        out <- data.frame(
          report_id = NA_character_,
          path = path_data[[outdir]]
        )
      } else {
        # NOTE: Calculating all possible combinations is not very efficient, but might be good enough.
        # Need a more robust way of naming files for a regex-based method to work reliably.
        ids <- list(
          sample_id = unique(sample_meta[[outdir]]$sample_id),
          reference_id = unique(reference_meta[[outdir]]$ref_id),
          report_id = unique(sample_meta[[outdir]]$report_group_ids)
        )
        combinations <- do.call(expand.grid, ids[id_types]) 
        combinations$combined <- apply(combinations, MARGIN = 1, paste0, collapse = sep)
        variable_name_part <- sub(basename(path_data[[outdir]]), pattern = regex, replacement = '\\1')
        variable_cols <- combinations[match(variable_name_part, combinations$combined), id_types, drop = FALSE]
        out <- cbind(variable_cols, path = path_data[[outdir]])
      }
      
      # If report ID is included and is NA, check if it can be inferred from the report input directory name
      is_report_input_dir <- grepl(outdir, pattern = 'report_group_data/.+_inputs$')
      if ("report_id" %in% id_types && is_report_input_dir) {
        out$report_id[is.na(out$report_id)] <- sub(outdir, pattern = '^.*report_group_data/(.+)_inputs$', replacement = '\\1')
      }
    }
    
    rownames(out) <- NULL
    tibble::as_tibble(out)
  })
  names(output) <- names(path_data)
  
  return(output)
}
