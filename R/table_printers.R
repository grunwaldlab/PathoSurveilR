#' Print a table of sample metadata
#'
#' Selects columns in the sample metadata to print and format the result for
#' use as a static table in PDF or an interactive table in HTML.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sample_meta_parsed()] output.
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @param ... Passed to `DT::datatable`.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sample_meta_table(path)
#' sample_meta_table(path, interactive = TRUE)
#'
#' @export
sample_meta_table <- function(input, interactive = FALSE, ...) {
  # Parse the input if it is a file/folder path
  if (inherits(input, 'data.frame')) {
    sample_data <- input
  } else {
    sample_data <- find_ps_data(input, target = 'sample_metadata', simplify = TRUE)
  }
  
  # Remove unneeded columns
  cols_to_remove <- c('outdir_path')
  sample_data <- sample_data[! colnames(sample_data) %in% cols_to_remove]

  # Print table
  if (interactive) {
    output <- DT::datatable(sample_data, class = "display nowrap", ...)
    output <- DT::formatStyle(output, colnames(sample_data), "white-space" = "nowrap")
  } else {
    output <- print_static_table(sample_data, compressed_cols = c('Forward Reads', 'Reverse Reads', 'Reference'))
  }

  return(output)
}


#' Get table of pipeline status data
#'
#' Return a formatted interactive table with the data on the issues encountered
#' by the pipeline, one row for each issue. The contents of all status message
#' files found in the given paths will be combined.
#'
#' @param input The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param summarize_by What variable to summarize results by, if any. Can be one
#'   of `'sample'`, `'message'`, `'workflow'`, or `NULL`. By default, all values
#'   in the message data is shown on its own row.
#' @param interactive Whether to produce interactive tables
#'   (TRUE) or static tables (FALSE). Defaults to TRUE if the environment
#'   supports HTML output, otherwise FALSE. Interactive tables offer enhanced
#'   browsing capabilities, while static tables are best for printed pdf
#'   reports.
#' @param ... Passed to `DT::datatable` for interactive output.
#'
#' @return A table with details for errors, warnings, notes.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' status_message_table(path)
#' status_message_table(path, interactive = TRUE)
#' status_message_table(path, summarize_by = 'sample')
#' status_message_table(path, summarize_by = 'workflow')
#' status_message_table(path, summarize_by = 'message')
#'
#' @export
status_message_table <- function(input, summarize_by = NULL, interactive = FALSE, ...) {
  # Check parameters
  valid_summarize_by <- c('sample', 'message', 'workflow')
  if (! is.null(summarize_by) && ! summarize_by %in% valid_summarize_by) {
    stop('The `summarize_by` parameter must be NULL or one of: ',
         paste0(valid_summarize_by, collapse = ', '))
  }

  # Find and parse message data
  message_data <- find_ps_data(input, target = 'messages', simplify = TRUE)

  # Make step status column
  symbol_key <- c(
    NOTE = '\u2705',
    WARNING = '\uD83D\uDFE1',
    ERROR = '\u274C'
  )
  symbol_key <- factor(symbol_key, levels = symbol_key, ordered = TRUE)
  message_data$status <- symbol_key[message_data$level]

  # Summarize by sample
  summarize_by_factor <- function(table, by) {
    if (nrow(table) == 0) {
      table <- data.frame(
        x = character(0),
        notes = numeric(0),
        warnings = numeric(0),
        errors = numeric(0),
        status = factor(0)
      )
    } else {
      table <- do.call(rbind, lapply(split(table, table[[by]]), function(part) {
        data.frame(
          x = unique(part[[by]]),
          notes = sum(part$level == 'NOTE'),
          warnings = sum(part$level == 'WARNING'),
          errors = sum(part$level == 'ERROR'),
          status = max(part$status)
        )
      }))
    }
    colnames(table)[1] <- by
    rownames(table) <- NULL
    return(table)
  }
  if (! is.null(summarize_by)) {
    if (summarize_by == 'sample') {
      message_data <- summarize_by_factor(message_data, 'sample_id')
    }
    if (summarize_by == 'workflow') {
      message_data <- summarize_by_factor(message_data, 'workflow')
    }
    if (summarize_by == 'message') {
      message_data <- summarize_by_factor(message_data, 'message')
      message_data$count <- rowSums(message_data[, c('notes', 'warnings', 'errors')])
      message_data <- message_data[, c('message', 'count', 'status')]
    }
  }

  # Format the analysis step to look nicer
  make_text_nice <- function(x) {
    x <- tools::toTitleCase(gsub(tolower(x), pattern = '_+', replacement = ' '))
    first_letter <- substr(x, 1, 1)
    vapply(seq_len(length(x)), FUN.VALUE = character(1), function(i) {
      sub(x[i], pattern = '^.?', replacement = toupper(first_letter[i]))
    })
  }
  columns_to_prettify <- c('workflow', 'level')
  columns_to_prettify <- columns_to_prettify[columns_to_prettify %in% colnames(message_data)]
  message_data[columns_to_prettify] <- lapply(message_data[columns_to_prettify], make_text_nice)

  # Sort by how bad the status is
  message_data <- message_data[order(message_data$status, decreasing = TRUE), ]

  # Reorder columns
  col_name_key <- c(
    status = 'Status',
    level = 'Type',
    sample_id = 'Sample ID',
    reference_id = 'Reference ID',
    workflow = 'Pipeline Step',
    count = 'Count',
    message = 'Message',
    notes = 'Notes',
    warnings = 'Warnings',
    errors = 'Errors'
  )
  col_name_key <- col_name_key[names(col_name_key) %in% colnames(message_data)]
  message_data <- message_data[, names(col_name_key)]

  # Rename columns
  print_data <- message_data
  colnames(print_data) <- col_name_key

  # Create pretty table to print
  if (interactive) {
    print_data$Status <- NULL
    output <- DT::datatable(print_data, options = list(pageLength = 5, autoWidth = TRUE),
                            escape = FALSE, rownames = as.character(message_data$status), ...)
    output <- DT::formatStyle(output, colnames(print_data), "white-space" = "nowrap")
  } else {
    rownames(print_data) <- NULL
    colnames(print_data)[colnames(print_data) == 'Status'] <- ''
    output <- print_static_table(print_data)
  }

  return(output)
}


#' Print best ANI match table
#'
#' Prints a table with the highest ANI matches for each sample given a pairwise
#' distance matrix.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param interactive Whether to produce interactive tables (TRUE) or static
#'   tables (FALSE). Defaults to TRUE if the environment supports HTML output,
#'   otherwise FALSE. Interactive tables offer enhanced browsing capabilities,
#'   while static tables are best for printed pdf reports.
#' @param ... Passed to `DT::datatable` for interactive output.
#' 
#' @export
estimated_ani_match_table <- function(path, interactive = FALSE, ...) {
  # Get best match table
  output <- make_best_match_table(
    pairwise_matrix = find_ps_data(path, target = 'ani_matrix_csv', simplify = TRUE),
    sample_data = find_ps_data(path, target = 'sample_metadata', simplify = TRUE),
    ref_data = find_ps_data(path, target = 'reference_metadata', simplify = TRUE)
  )

  # Print table
  output$best_ref_value <- ifelse(is.na(output$best_ref_value),
                                          "NA", format_number(output$best_ref_value))
  output$best_sample_value <- ifelse(is.na(output$best_sample_value),
                                             "NA", format_number(output$best_sample_value))
  colnames(output) <- c('Sample', 'Closest reference', 'Reference ANI (%)', 'Closest sample', 'Sample ANI (%)')
  if (interactive) {
    output <- DT::datatable(output, class = "display nowrap", ...)
    output <- DT::formatStyle(output, colnames(output), "white-space" = "nowrap")
  } else {
    output <- print_static_table(output)
  }
  return(output)
}


#' Print best POCP table
#'
#' Prints a table with the highest POCP for each sample given a pairwise
#' distance matrix.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param interactive Whether to produce interactive tables (TRUE) or static
#'   tables (FALSE). Defaults to TRUE if the environment supports HTML output,
#'   otherwise FALSE. Interactive tables offer enhanced browsing capabilities,
#'   while static tables are best for printed pdf reports.
#' @param ... Passed to `DT::datatable` for interactive output.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_match_table(path)
#' pocp_match_table(path, interactive = TRUE)
#'
#' @export
pocp_match_table <- function(path, interactive = FALSE, ...) {
  # Get best match table
  output <- make_best_match_table(
    pairwise_matrix = find_ps_data(path, target = 'ani_matrix_csv', simplify = TRUE),
    sample_data = find_ps_data(path, target = 'sample_metadata', simplify = TRUE),
    ref_data = find_ps_data(path, target = 'reference_metadata', simplify = TRUE)
  )

  # Print table
  output$best_ref_value <- ifelse(is.na(output$best_ref_value),
                                          "NA", format_number(output$best_ref_value))
  output$best_sample_value <- ifelse(is.na(output$best_sample_value),
                                             "NA", format_number(output$best_sample_value))
  colnames(output) <- c('Sample', 'Closest reference', 'Reference POCP (%)', 'Closest sample', 'Sample POCP (%)')
  if (interactive) {
    output <- DT::datatable(output, class = "display nowrap", ...)
    output <- DT::formatStyle(output, colnames(output), "white-space" = "nowrap")
  } else {
    output <- print_static_table(output)
  }
  return(output)
}


#' Print best match table
#'
#' Prints a table with the highest value for each sample given a pairwise
#' distance matrix.
#'
#' @param pairwise_matrices A `list` of `matrix` with all pairwise comparisons
#'   between samples and references
#' @param sample_data The sample metadata
#' @param ref_data The reference metadata
#'
#' @keywords internal
make_best_match_table <- function(pairwise_matrix, sample_data, ref_data) {
  sample_ids <- sample_data$sample_id[sample_data$sample_id %in% colnames(pairwise_matrix)]
  ref_ids <- ref_data$ref_id[ref_data$ref_id %in% colnames(pairwise_matrix)]
  do.call(rbind, lapply(sample_ids, function(id) { # loop over sample IDs and combine results into a table
    best_ref_match_name <- NA
    best_ref_match_value <- NA
    if (length(ref_ids) >= 1) {
      ref_samp_comp <- pairwise_matrix[id, colnames(pairwise_matrix) %in% ref_data$ref_id, drop = FALSE]
      if (! all(is.na(ref_samp_comp))) {
        best_ref_match_id <- colnames(ref_samp_comp)[which.max(ref_samp_comp)]
        best_ref_match_value <- unname(unlist(ref_samp_comp[id, best_ref_match_id]))
        best_ref_match_name <- ref_data$ref_name[ref_data$ref_id == best_ref_match_id]
      }
    } 
    best_samp_match_name <- NA
    best_samp_match_value <- NA
    if (length(sample_ids) > 1) {
      samp_samp_comp <- pairwise_matrix[id, colnames(pairwise_matrix) %in% sample_ids & colnames(pairwise_matrix) != id, drop = FALSE]
      if (! all(is.na(samp_samp_comp))) {
        best_samp_match_id <- colnames(samp_samp_comp)[which.max(samp_samp_comp)]
        best_samp_match_value <- unname(unlist(samp_samp_comp[id, best_samp_match_id]))
        best_samp_match_name <- sample_data$name[sample_data$sample_id == best_samp_match_id]
      }
    } 
    data.frame(
      check.names = FALSE,
      sample_name = sample_data$name[sample_data$sample_id == id],
      best_ref = best_ref_match_name,
      best_ref_value = best_ref_match_value,
      best_sample = best_samp_match_name,
      best_sample_value = best_samp_match_value
    )
  }))
}


#' Create an summary table from sendsketch output table
#'
#' This function processes sendsketch data from pathogen surveillance nextflow
#' pipeline and generates an interactive HTML table displaying the top hits from
#' the sketch-based analysis. The resulting table includes metrics such as
#' Weighted Kmer IDentity (WKID), Average Nucleotide Identity (ANI), and
#' Completeness, providing insights into genomic similarity and representational
#' completeness between the query and reference genomes. The table visually
#' encodes the percentage values using horizontal bars for an intuitive and
#' accessible presentation.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sendsketch_parsed()] output. Only the best hits are returned, based on the
#'   default behavior of [sendsketch_best_hits()]. To change this behavior, pass
#'   in the results of running [sendsketch_best_hits()].
#' @param interactive Whether or not to produce an interactive
#'   HTML/javascript-based table or a static one.
#'
#' @return An interactive DT::datatable object that renders as an HTML table
#'   when used in an R Markdown document or Shiny application. The table will
#'   have interactive features such as sorting and search enabled.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_table(path)
#' sendsketch_table(path, interactive = TRUE)
#'
#' @export
sendsketch_table <- function(path, interactive = FALSE) {
  # Parse the input if it is a file/folder path
  if (is.data.frame(path)) {
    sketch_data <- path
  } else {
    sketch_data <- sendsketch_parsed(path, only_best = TRUE)
  }

  # Sort and filter data
  final_table <- sketch_data[, c('sample_id', 'WKID', 'ANI', 'Complt', 'taxName'), drop = FALSE]
  new_col_names <- c('Sample', 'WKID (%)', 'ANI (%)', 'Completeness (%)', 'Top Hit')

  if (interactive) {
    # Define a function called 'bordered_bar'
    bordered_bar <- function(value, color) {
      sprintf('<div style="border: 1px solid gray; width: 100%%; border-radius: 12px;">
              <div style="width: %s%%; background-color: %s; border-radius: 12px; text-align: center;">%s</div>
             </div>',
              value, color, value)
    }

    # Apply the bordered_bar function
    final_table$WKID <- vapply(final_table$WKID, FUN.VALUE = character(1), function(x) bordered_bar(x, 'lightblue'))
    final_table$ANI <- vapply(final_table$ANI, FUN.VALUE = character(1), function(x) bordered_bar(x, 'lightgreen'))
    final_table$Complt <- vapply(final_table$Complt, FUN.VALUE = character(1), function(x) bordered_bar(x, 'lightpink'))

    # Render the table using DT for HTML output
    names(final_table) <- new_col_names
    return(
      DT::datatable(final_table,
                    options = list(
                      pageLength = 10,
                      autoWidth = TRUE,
                      columnDefs = list(
                        list(width = '150px', targets = c(1, 2, 3))
                      ),
                      searchHighlight = TRUE
                    ),
                    rownames = FALSE,
                    escape = FALSE
      )
    )
  } else {
    names(final_table) <- new_col_names
    return(print_static_table(final_table))
  }
}