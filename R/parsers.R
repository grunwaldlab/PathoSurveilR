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
      sketch_data <- make_empty_data_frame(
        c("sample_id", "WKID", "KID", "ANI", "SSU", "SSULen", "Complt", 
          "Contam", "Contam2", "uContam", "Score", "E-Val", "Depth", "Depth2", 
          "Volume", "RefHits", "Matches", "Unique", "Unique2", "Unique3", 
          "noHit", "Length", "TaxID", "ImgID", "gBases", "gKmers", "gSize", 
          "gSeqs", "GC", "rDiv", "qDiv", "rSize", "qSize", "cHits", "taxName", 
          "file", "seqName", "taxonomy", "outdir_path")
      )
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
      class_data <- lookup_ncbi_taxon_id_classification(sketch_data$TaxID, match_input = FALSE, simplify = FALSE)
      classifications  <- vapply(class_data, FUN.VALUE = character(1), function(x) {
        paste0(x$rank, ':', x$name, collapse = ';')
      })
      sketch_data$taxonomy <- classifications[sketch_data$TaxID]
    }
    
    tibble::as_tibble(sketch_data)
  })
  
  postprocess_table_list(output, simplify = simplify)
}


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
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param ranks_as_cols If `TRUE`, return a table with a column for each rank
#'   instead of a table with one taxon per row.
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
sendsketch_taxonomy_parsed <- function(path, as_table = TRUE, only_shared = FALSE,
                                       sample_id = NULL, only_best = FALSE, simplify = TRUE) {
  sendsketch_data <- sendsketch_parsed(path, sample_id = sample_id, only_best = only_best,
                                       update_taxonomy = FALSE, simplify = FALSE)
  all_tax_ids <- unlist(lapply(sendsketch_data, function(x) x$TaxID))
  class_data <- lookup_ncbi_taxon_id_classification(all_tax_ids, match_input = FALSE, simplify = FALSE)
  
  process_one <- function(table) {
    if (nrow(table) == 0) {
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
    
    table_tax_ids <- unique(table$TaxID)
    table_class_data <- class_data[table_tax_ids]
    
    # Remove data for ranks not shared among all results
    if (only_shared) {
      if (simplify) {
        all_ranks <- unlist(lapply(class_data, function(x) unique(x$rank)))
        is_shared <- table(all_ranks) == length(class_data)
      } else {
        all_ranks <- unlist(lapply(table_class_data, function(x) unique(x$rank)))
        is_shared <- table(all_ranks) == length(table_class_data)
      }
      shared_ranks <- names(is_shared)[is_shared]
      table_class_data <- lapply(table_class_data, function(x) x[x$rank %in% shared_ranks, , drop = FALSE])
    }
    
    # Remove ranks that are placeholders for missing ranks
    placeholder_ranks <- c('clade', 'no rank')
    table_class_data <- lapply(table_class_data, function(x) x[! x$rank %in% placeholder_ranks, , drop = FALSE])
    
    if (as_table) {
      # Reformat with ranks as columns
      parts <- lapply(table_tax_ids, function(id) {
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
  
  postprocess_table_list(lapply(sendsketch_data, process_one), simplify = simplify)
}


#' Get taxa predicted by sendsketch results
#'
#' Return which taxa are predicted to be present given sendsketch results
#' associated with directory paths containing `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param only_best Only return taxa predicted by the best hit for each
#'   combination of report group and sample.
#' @inheritParams sendsketch_parsed
#' @inheritParams postprocess_table_list
#'
#' @return A [tibble::tibble()] with one row per taxon found
#' @family parsers
#'
#' @export
sendsketch_taxa_present <- function(path, ani_thresh = c(species = 95, genus = 90, family = 70),
                                    comp_thresh = c(species = 40, genus = 15, family = 5),
                                    sample_id = NULL, only_best = FALSE, simplify = TRUE) {
  all_tax_data <- sendsketch_taxonomy_parsed(path = path, as_table = FALSE, only_shared = FALSE,
                                         sample_id = sample_id, only_best = only_best, simplify = FALSE)
  all_sketch_data <- sendsketch_parsed(path = path, sample_id = sample_id, only_best = only_best, simplify = FALSE)
  
  output <- lapply(names(all_tax_data), function(out_path) {
    tax_data <- all_tax_data[[out_path]]
    sketch_data <- all_sketch_data[[out_path]]
    filter_and_extract <- function(rank) {
      subset_ids <- sketch_data$TaxID[sketch_data$ANI > ani_thresh[rank] & sketch_data$Complt > comp_thresh[rank]]
      subset_class_data <- tax_data[tax_data$query_id %in% subset_ids & tax_data$rank == rank, , drop = FALSE]
      subset_class_data[, c('id', 'name', 'rank'), drop = FALSE]
    }
    output_data <- unique(do.call(rbind, lapply(names(ani_thresh), filter_and_extract)))
  })
  names(output) <- names(all_tax_data)

  postprocess_table_list(output, simplify = simplify)
}


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
#' @param taxon_id The NCBI taxon IDs of families to return data for.
#' @inheritParams postprocess_table_list
#'
#' @return A table of reference metadata
#' 
#' @export
considered_ref_meta_parsed <- function(path = NULL, taxon_id = NULL, simplify = TRUE) {
  # Check input type parameters
  if (! is.null(taxon_id) && ! is.character(taxon_id)) {
    stop('The argument `taxon_id` must be a character vector')
  }
  
  # Find data on input files and filter by 'taxon_id' argument
  all_path_data <- considered_ref_meta_path_data(path, simplify = FALSE)
  if (! is.null(taxon_id)) {
    all_path_data <- lapply(all_path_data, function(x) x[x$taxon_id %in% taxon_id, , drop = FALSE])
  }
  
  output <- lapply(all_path_data, function(path_data) {
    assem_data <- do.call(rbind, lapply(path_data$path, function(path) {
      out <- read.table(path, header = TRUE, sep = '\t', comment.char = '')
      family_id <- gsub(basename(path), pattern = '.tsv', replacement = '', fixed = TRUE)
      if (nrow(out) > 0) {
        out$family_id <- family_id
      }
      return(out)
    }))
    
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
  })
  
  postprocess_table_list(output, simplify = simplify)
}
