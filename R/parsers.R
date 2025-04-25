#' Get parsed pipeline status data
#'
#' Return a [tibble::tibble()] (table) with the status messages produced by the
#' pathogensuriveillance pipeline. The contents of all status message files
#' found in the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
#'
#' @return A [tibble::tibble()] with the messages from all input paths if
#'   `simplify = TRUE` or a list of such [tibble::tibble()]s for each output
#'   directory found if `simplify = FALSE`.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_parsed(path)
#'
#' @export
status_message_parsed <- function(path, simplify = TRUE) {
  path_data <- status_message_path_data(path, simplify = FALSE)
  output <- lapply(path_data, function(table) {
    if (nrow(table) > 0) {
      do.call(rbind, lapply(table$path, utils::read.table, sep = '\t', check.names = FALSE, header = TRUE))
    } else {
      tibble::tibble(
        report_id = character(0),
        sample_id = character(0),
        reference_id = character(0),
        workflow = character(0),
        level = character(0),
        message = character(0)
      )
    }
  })
  postprocess_table_list(output, simplify = simplify)
}


#' Get summary of pipeline status data
#'
#' Return a [tibble::tibble()] (table) with the numbers of issues encountered by the
#' pipeline. The contents of all status message files found in the given paths
#' will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
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
status_message_parsed_summary <- function(path, simplify = TRUE) {
  message_data <- status_message_parsed(path, simplify = FALSE)
  output <- lapply(message_data, function(table) {
    if (nrow(table) == 0) {
      tibble::tibble(
        workflow = character(0),
        errors =  numeric(0),
        warnings = numeric(0),
        notes = numeric(0),
        samples = numeric(0),
        references = numeric(0),
        report_groups = numeric(0)
      )
    } else {
      do.call(rbind, lapply(split(table, table$workflow), function(subset) {
        tibble::tibble(
          workflow = unique(subset$workflow),
          errors = sum(subset$level == 'ERROR'),
          warnings = sum(subset$level == 'WARNING'),
          notes = sum(subset$level == 'NOTE'),
          samples = length(unique(subset$sample_id)),
          references = length(unique(subset$reference_id)),
          report_groups = length(unique(subset$report_group_id))
        )
      }))
    }
  })
  postprocess_table_list(output, simplify = simplify)
}


#' Get parsed sample metadata
#'
#' Return a [tibble::tibble()] (table) with sample metadata. The contents of all
#' sample metadata files found in the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
#'
#' @return A [tibble::tibble()] with the sample metadata if
#'   `simplify = TRUE` or a list of such [tibble::tibble()]s for each output
#'   directory found if `simplify = FALSE`.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_parsed(path)
#'
#' @export
sample_meta_parsed <- function(path, simplify = TRUE) {
  meta_paths <- sample_meta_path(path, simplify = FALSE)
  output <- lapply(meta_paths, function(p) {
    x <- do.call(combine_data_frames, lapply(p, utils::read.csv, check.names = FALSE, sep = '\t'))
    x[] <- lapply(x, function(col_data) ifelse(col_data == 'null', NA_character_, col_data))
    tibble::as_tibble(x)
  })
  postprocess_table_list(output, simplify = simplify)
}


#' Get parsed reference metadata
#'
#' Return a [tibble::tibble()] (table) with reference metadata. The contents of
#' all reference metadata files found in the given paths will be combined. This
#' only contains data about references used by at least one step of the pipeline.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
#'
#' @return A [tibble::tibble()] with the reference metadata if
#'   `simplify = TRUE` or a list of such [tibble::tibble()]s for each output
#'   directory found if `simplify = FALSE`.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' ref_meta_parsed(path)
#'
#' @export
ref_meta_parsed <- function(path, simplify = TRUE) {
  meta_paths <- ref_meta_path(path, simplify = FALSE)
  output <- lapply(meta_paths, function(p) {
    x <- do.call(combine_data_frames, lapply(p, utils::read.csv, check.names = FALSE, sep = '\t'))
    x[] <- lapply(x, function(col_data) ifelse(col_data == 'null', NA_character_, col_data))
    tibble::as_tibble(x)
  })
  postprocess_table_list(output, simplify = simplify)
}


#' Parse distance matrices
#'
#' Parse distance matrices
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param path_func A path finder function
#' @param simplify If `TRUE`, combine output into a single matrix, filling in
#'   the gaps with `NA`.
#' @param ... passed to read.table
#' @inheritParams combine_matrices
#' 
#' @return a `list` of matrices
#'
#' @keywords internal
generic_matrix_parsed <- function(path, path_func, simplify = TRUE, as_edge_list = FALSE, ...) {
  matrix_paths <- path_func(path, simplify = FALSE)
  output <- lapply(matrix_paths, function(paths) {
    matrices <- lapply(paths, function(one_path) {
      mat <- as.matrix(utils::read.table(one_path, check.names = FALSE, header = TRUE, ...))
      rownames(mat) <- colnames(mat)
      mat[mat == 0] <- NA
      return(mat)
    })
    combine_matrices(matrices, as_edge_list = as_edge_list)
  })
  names(output) <- matrix_paths
  
  if (simplify) {
    if (as_edge_list) {
      output <- do.call(rbind, output)
      output <- unique(output)
      rownames(output) <- NULL
    } else {
      output <- combine_matrices(output, as_edge_list = FALSE)
    }
  }
  
  return(output)
}


#' Get estimated ANI distance matrix
#'
#' Return a list of [base::matrix()] with the estimated pairwise ANI values
#' calculated by sourmash.
#'
#' @inheritParams generic_matrix_parsed
#' 
#' @return a `list` of ANI matrices or edge lists depending on the values of `simplify` and `as_edge_list`
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_matrix_parsed(path)
#'
#' @export
estimated_ani_matrix_parsed <- function(path, simplify = FALSE, as_edge_list = FALSE) {
  generic_matrix_parsed(path, estimated_ani_matrix_path, simplify = simplify, as_edge_list = as_edge_list, sep = ',')
}


#' Get POCP matrix
#'
#' Return a list of [base::matrix()] with the POCP values based on a core
#' gene analysis using pirate.
#'
#' @inheritParams generic_matrix_parsed
#' 
#' @return a `list` of ANI matrices or edge lists depending on the values of `simplify` and `as_edge_list`
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_matrix_parsed(path)
#'
#' @export
pocp_matrix_parsed <- function(path, simplify = FALSE, as_edge_list = FALSE) {
  generic_matrix_parsed(path, pocp_matrix_path, simplify = simplify, as_edge_list = as_edge_list, sep = '\t')
}


#' Get parsed report groups
#'
#' Return a [tibble::tibble()] (table) with report group ID. The groups found in
#' the given paths will be combined.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
#'
#' @return A [base::character()] vector with the groups from all input paths
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' run_info_parsed(path)
#'
#' @export
run_info_parsed <- function(path, simplify = TRUE) {
  path_data <- run_info_path_data(path, simplify = FALSE)
  output <- lapply(path_data, function(table) {
    if (nrow(table) == 0) {
      tibble::tibble(
        command_line = character(0),
        commit_id = character(0),
        container_engine = character(0),
        profile = character(0),
        revision = character(0),
        run_name = character(0),
        session_id = character(0),
        start_time = character(0),
        nextflow_version = character(0),
        pipeline_version = character(0)
      )
    } else {
      do.call(rbind, lapply(table$path, function(p) {
        run_info <- yaml::read_yaml(p, readLines.warn = FALSE)
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
  })
  postprocess_table_list(output, simplify = simplify)
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
#' @param update_taxonomy If `TRUE` look up the current taxonomy classification
#'   using the NCBI taxon ID rather than using the existing one.
#' @inheritParams sendsketch_path_data
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
sendsketch_parsed <- function(path, sample_id = NULL, only_best = FALSE, update_taxonomy = FALSE, simplify = TRUE) {
  path_data <- sendsketch_path_data(path, simplify = FALSE)
  
  output <- lapply(path_data, function(table) {
    if (! is.null(sample_id)) {
      table <- table[table$sample_id %in% sample_id]    
    }
    
    # Parse all files and combine them into one data frame
    if (nrow(table) == 0) {
      cols <- c("sample_id", "WKID", "KID", "ANI", "SSU", "SSULen", "Complt", 
        "Contam", "Contam2", "uContam", "Score", "E-Val", "Depth", "Depth2", 
        "Volume", "RefHits", "Matches", "Unique", "Unique2", "Unique3", 
        "noHit", "Length", "TaxID", "ImgID", "gBases", "gKmers", "gSize", 
        "gSeqs", "GC", "rDiv", "qDiv", "rSize", "qSize", "cHits", "taxName", 
        "file", "seqName", "taxonomy", "outdir_path")
      sketch_data <-  data.frame(matrix(vector(), 0, length(cols), dimnames=list(c(), cols)))
    } else {
      sketch_data <- do.call(rbind, lapply(seq_len(nrow(table)), function(i) {
        data <- suppressWarnings(utils::read.table(table$path[i], sep = '\t', skip = 2, check.names = FALSE, header = TRUE))
        return(cbind(
          sample_id = rep(table$sample_id[i], nrow(data)),
          data,
          outdir_path = table$outdir_path[i]
        ))
      }))
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
      all_tax_ids <- unique(sketch_data$TaxID)
      tax_id_batches <- split(all_tax_ids, ceiling(seq_along(all_tax_ids) / 50))
      classifications <- unlist(lapply(tax_id_batches, function(tax_ids) {
        raw_xml = rentrez::entrez_fetch(db = 'taxonomy', id = tax_ids, rettype = 'xml')
        Sys.sleep(0.33)
        get_capture_group <- function(x, pattern) {
          matches <- regmatches(x, gregexpr(pattern, x))[[1]]
          sub(matches, pattern = pattern, replacement = '\\1')
        }
        class_xml <- get_capture_group(raw_xml, '<LineageEx>(.+?)</LineageEx>')
        class_names <- lapply(class_xml, get_capture_group, pattern = '<ScientificName>(.+?)</ScientificName>')
        class_ranks <- lapply(class_xml, get_capture_group, pattern = '<Rank>(.+?)</Rank>')
        classifications  <- vapply(seq_len(length(class_names)), FUN.VALUE = character(1), function(i) {
          paste0(class_ranks[[i]], ':', class_names[[i]], collapse = ';')
        })
      }))
      names(classifications) <- all_tax_ids
      sketch_data$taxonomy <- classifications[sketch_data$TaxID]
    }
    
    tibble::as_tibble(sketch_data)
  })
  
  postprocess_table_list(output, simplify = simplify)
}


#' Get parsed sendsketch taxonomy
#'
#' Return the taxonomic classification in sendsketch results associated with
#' directory paths containing `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param as_table If `TRUE`, return a table of classification data with a
#'   column for each rank rather than a character vector with all of the
#'   classifications.
#' @param remove_ranks If `TRUE`, remove the rank information from the taxonomy.
#'   Only has an effect if `as_table` is `FALSE`.
#' @param only_shared If `TRUE`, only return the ranks that are present in all
#'   of the inputs. Only has an effect if `as_table` is `TRUE`.
#' @param append If `TRUE`, append output to the rest of the parsed sendskech
#'   results. Only has an effect if `as_table` is `TRUE`.
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
sendsketch_taxonomy_parsed <- function(path, as_table = TRUE, remove_ranks = FALSE, only_shared = FALSE, append = FALSE, 
                                       sample_id = NULL, only_best = FALSE, update_taxonomy = FALSE, simplify = TRUE) {
  sendsketch_data <- sendsketch_parsed(path, sample_id = sample_id, only_best = only_best, update_taxonomy = update_taxonomy, simplify = FALSE)
  
  if (as_table) {
    output <- lapply(sendsketch_data, function(table) {
      # Reformat with ranks as columns
      parts <- lapply(seq_len(nrow(table)), function(i) {
        split_tax <- strsplit(table$taxonomy[i], split = ';', fixed = TRUE)[[1]]
        taxon_names <- gsub(split_tax, pattern = '^[a-z]+:', replacement = '')
        ranks <- gsub(split_tax, pattern = '^([a-z]*):?.+$', replacement = '\\1')
        names(taxon_names) <- ifelse(ranks == '', 'tip', ranks)
        tibble::as_tibble(as.list(c(sample_id = names(table$taxonomy)[i], taxon_names)))
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
      
      # Remove columns not shared by all samples
      if (only_shared) {
        col_has_na <- apply(output, MARGIN = 2, function(col) any(is.na(col)))
        output <- output[, ! col_has_na]
      }
      
      # Add back on to input data
      if (append) {
        output <- cbind(table, output)
      }
      
      tibble::as_tibble(output)
    })
    output <- postprocess_table_list(output, simplify = simplify)
  } else {
    output <- lapply(sendsketch_data, function(table) {
      classifications <- table$taxonomy
      if (remove_ranks) {
        classifications <- gsub(classifications, pattern = ';?[a-z]+:', replacement = ';')
        classifications <- gsub(classifications, pattern = '^;', replacement = '')
      }
      names(classifications) <- table$sample_id
      classifications[! duplicated(paste0(classifications, names(classifications)))]
    })
    if (simplify) {
      out_names <- unlist(lapply(output, names), use.names = FALSE)
      output <- unlist(output, use.names = FALSE)
      names(output) <- out_names
    }
  }

  return(output)
}


# #' Get taxa predicted by sendsketch results
# #'
# #' Return which taxa are predicted to be present given sendsketch results
# #' associated with directory paths containing `pathogensurveillance` output.
# #'
# #' @param path The path to one or more folders that contain pathogensurveillance
# #'   output.
# #' @param only_best Only return taxa predicted by the best hit for each
# #'   combination of report group and sample.
# #' @inheritParams sendsketch_parsed
# #' @inheritParams postprocess_table_list
# #'
# #' @return A [tibble::tibble()] with one row per taxon found
# #' @family parsers
# #'
# #' @export
# sendsketch_taxa_present <- function(path, ani_thresh = c(species = 95, genus = 90, family = 70),
#                                     comp_thresh = c(species = 40, genus = 15, family = 5),
#                                     sample_id = NULL, only_best = FALSE, simplify = TRUE) {
#   sketch_data <- sendsketch_taxonomy_parsed(path, sample_id = sample_id, append = TRUE, simplify = FALSE)
#   
#   lapply(sketch_data, function(table) {
#     # Filter data by threshold and extract passing taxon names
#     filter_and_extract <- function(rank) {
#       subset_ids <- table$TaxID[table$ANI > ani_threshold[rank] & table$Complt > complt_threshold[rank]]
#       subset_class_data <- class_data[class_data$tip_taxon_id %in% subset_ids & class_data$rank == rank, , drop = FALSE]
#       subset_class_data[, c('taxon_id', 'name', 'rank'), drop = FALSE]
#     }
#     output_data <- unique(do.call(rbind, lapply(names(ani_threshold), filter_and_extract)))
#     
#   })
#   
# }


#' Parse Software Version Metadata from YAML File
#'
#' Reads a YAML file containing software version information and transforms it
#' into a tibble. Each software module and program, along with its version, is
#' extracted and organized.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
#' @return A [tibble::tibble()] with columns `module`, `program`, and `version`
#'   detailing the software versions.
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' software_version_parsed(path)
#'
#' @export
software_version_parsed <- function(path, simplify = TRUE) {
  version_paths <- software_version_path(path, simplify = FALSE)
  parse_one <- function(version_path) {
    raw_version_data <- yaml::read_yaml(version_path)
    version_data <- tibble::tibble(
      module = rep(names(raw_version_data), vapply(raw_version_data, FUN.VALUE = numeric(1), length)),
      program = unname(unlist(lapply(raw_version_data, names))),
      version = unname(unlist(raw_version_data))
    )
    return(version_data)
  }
  output <- lapply(version_paths, function(p) {
    unique(do.call(rbind, lapply(p, parse_one)))
  })
  output <- postprocess_table_list(output, simplify = simplify)
  if (simplify) {
    output <- unique(output)
  }
  return(output)
}


#' Get core gene phylogeny
#'
#' Return a list of [ape::phylo()] objects named by file path by searching for
#' folders containing pathogensurveillance output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @inheritParams postprocess_table_list
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' core_tree_parsed(path)
#'
#' @export
core_tree_parsed <- function(path, simplify = TRUE) {
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
#' @inheritParams postprocess_table_list
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_tree_parsed(path)
#'
#' @export
variant_tree_parsed <- function(path, rename = TRUE, simplify = TRUE) {
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
#' @inheritParams postprocess_table_list
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' busco_tree_parsed(path)
#'
#' @export
busco_tree_parsed <- function(path, simplify = TRUE) {
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
#' @inheritParams postprocess_table_list
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' multigene_tree_parsed(path)
#'
#' @export
multigene_tree_parsed <- function(path, simplify = TRUE) {
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
#' @inheritParams postprocess_table_list
#'
#' @return List of [ape::phylo()] objects named by file path
#' @family parsers
#'
#' @keywords internal
tree_parsed <- function(path, simplify = TRUE) {
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
#' @inheritParams postprocess_table_list
#' 
#' @return list of [ape::DNAbin()], named by alignment file
#' @family parsers
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_align_parsed(path)
#'
#' @export
variant_align_parsed <- function(path, rename = TRUE, simplify = TRUE) {
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
#' @inheritParams postprocess_table_list
#'
#' @return A table of reference metadata
#' 
#' @export
considered_ref_meta_parsed <- function(path = NULL, family = NULL, json_path = NULL, simplify = TRUE) {
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
        reference_id = gsub(assem_data$accession, pattern = '[\\/:*?"<>| .]', replacement = '_'),
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
      reference_id = character(0),
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
