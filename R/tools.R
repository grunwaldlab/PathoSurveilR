# This file contains internal functions that are generally useful in multiple contexts


#' Wrapper function to print static tables
#'
#' This wrapper function is meant to allow the method of printing tables used by
#' many functions to be changed in the future. This is primarily used for PDF
#' output.
#'
#' @param data The table to print.
#' @param compressed_cols The named of columns to attempt to make shorter by
#'   moving shared text at the beginning and end to footnotes and replacing with
#'   `...`
#' @param max_nchar The number of characters that the longest entry in a column
#'   can be before the `compressed_cols` option takes effect.
#'
#' @keywords internal
print_static_table <- function(data, compressed_cols = NULL, max_nchar = 20) {
  # Compress column contents
  if (!is.null(compressed_cols)) {
    # Dont compress column that already have small values
    compress_needed <- unlist(lapply(compressed_cols, function(col) {
      ! all(is.na(data[[col]])) && max(nchar(data[[col]])) > max_nchar
    }))
    compressed_cols <- compressed_cols[compress_needed]
    if (length(compressed_cols) > 0 ) {
      # Remove starts
      starts <- lapply(data[compressed_cols], shared_char)
      data[compressed_cols] <- lapply(compressed_cols, function(col_name) {
        column <- data[[col_name]]
        start <- starts[[col_name]]
        if (start != "") {
          column <- sub(column, pattern = paste0("^", start), replacement = '\u2026')
        }
        return(column)
      })
      # Remove ends
      ends <- lapply(data[compressed_cols], shared_char, end = TRUE)
      data[compressed_cols] <- lapply(compressed_cols, function(col_name) {
        column <- data[[col_name]]
        end <- ends[[col_name]]
        if (end != "") {
          column <- sub(column, pattern = paste0(end, "$"), replacement = '\u2026')
        }
        return(column)
      })
      # Create footnotes
      footnotes <- unlist(lapply(compressed_cols, function(col_name) {
        start <- starts[[col_name]]
        end <- ends[[col_name]]
        if (start != "" && end != "") {
          note <- paste0('All values in column "', col_name, '" start with "', start, '" and end with "', end, '".')
        } else if (start != "" ) {
          note <- paste0('All values in column "', col_name, '" start with "', start, '".')
        } else if (end != "" ) {
          note <- paste0('All values in column "', col_name, '" end with "', end, '".')
        } else {
          note <- NA_character_
        }
      }))
      names(footnotes) <- compressed_cols
      footnotes <- footnotes[! is.na(footnotes)]
      # Modify column names ( cant figure out how to get superscripts to render correctly )
      colnames(data)[colnames(data) %in% names(footnotes)] <- paste0(colnames(data)[colnames(data) %in% names(footnotes)], ' (', seq_along(footnotes), ')')
    } else {
      footnotes <- character(0)
    }
  } else {
    footnotes <- character(0)
  }

  # Print table
  is_multi_page_table <- nrow(data) > 45
  table <- kableExtra::kbl(data, booktabs = TRUE, longtable = is_multi_page_table) %>%
    kableExtra::kable_styling(full_width = FALSE, latex_options = c("hold_position", "repeat_header", "scale_down"))
  if (length(footnotes) > 0) {
    table <- kableExtra::footnote(table, number = footnotes)
  }
  return(table)
}


#' Identify conserved start/ends of characters
#'
#' Find any starts or ends of strings that are the same in all elements of a character vector
#'
#' @param col The character vector
#' @param end If `TRUE` then return shared end instead of shared start
#'
#' @keywords internal
shared_char <- function(col, end = FALSE) {

  reverse_char <- function(x) {
    split <- strsplit(x, split = "")
    reversed <- lapply(split, rev)
    return(unlist(lapply(reversed, paste, collapse = "")))
  }

  if (all(is.na(col))) {
    return("")
  }
  shorest_length <- min(unlist(lapply(col, nchar)), na.rm = TRUE)
  result <- ""
  if (end) {
    col <- reverse_char(col)
  }
  indexes <- 1:shorest_length
  for (index in indexes) {
    unique_starts <- unique(substr(col, start = 1, stop = index))
    if (length(unique_starts) == 1) {
      result <- unique_starts
    } else {
      break
    }
  }
  if (end) {
    result <- reverse_char(result)
  }
  return(result)
}



#' Format a number for printing
#'
#' Shortens a number to a specified number of significant digits.
#' Taken from https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
#'
#' @keywords internal
format_number <- function(nums, sig_fig = 4) {
  formatC(signif(nums, digits = sig_fig), digits = sig_fig, format="fg", flag="#")
}


#' Print figures with dropdown selector
#'
#' Prints HTML code to show a plot based on the value in a dropdown selector.
#' For use with Quarto/Rmarkdown/Knitr, put in a chunk with the option
#' `results='asis'`.
#'
#' @param plot_func A function to produce a plot using a single input from the
#'   `selector` input.
#' @param selector A named list of character vectors, that will be used to
#'   generate the input to `plot_func` to make plots. Plots will be made for
#'   every combination of values in the input. The number of character vectors
#'   should equal the number of arguments taken by `plot_func`.
#' @param id_prefix The prefix added to element IDs to distinguish this plot
#'   from others. This must be unique amoung other calls to this functions in a
#'   single HTML file.
#' @param imglist_class The CSS class used as the prefix for the IDs of the
#'   selector dropdown HTML elements. This must be unique amoung other calls to
#'   this functions in a single HTML file.
#' @param hide_single_selector If `TRUE`, don't show selector dropdown if it
#'   contains only one choice.
#' @param zoom If `TRUE`, add JS code to allow the image to be zoomed.
#' @param zoom_slider If `TRUE`, add a slider to the plot to control the zoom.
#' @param ... Passed to [grDevices::png()]
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' n <- c(10, 100, 1000)
#' base_plot_func <- function(x) {
#'   hist(rnorm(x))
#' }
#' print_figures_with_selector(base_plot_func, list(Count = n), 'base_test_id')
#'
#' ggplot_func <- function(x) {
#'   df <- data.frame(
#'     sex=factor(rep(c("F", "M"), each=x)),
#'     weight=round(c(rnorm(x, mean=55, sd=5), rnorm(x, mean=65, sd=5)))
#'   )
#'   print(ggplot(df, aes(x=weight)) + geom_histogram())
#' }
#' print_figures_with_selector(ggplot_func, list(Number = n), 'ggplot_test_id')
print_figures_with_selector <- function(plot_func, selector, id_prefix, imglist_class = paste0(id_prefix, '_list'), hide_single_selector = TRUE, zoom = TRUE, zoom_slider = FALSE, ...) {

  # Get combinations of input parameters
  plot_data <- expand.grid(selector, stringsAsFactors = FALSE)

  # Make plots encoded in base64
  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }
  plot_data$base64_plot <- unlist(lapply(seq_len(nrow(plot_data)), function(i) {
    temp_path <- tempfile(fileext = '.png')
    grDevices::png(temp_path, ...)
    args <- unname(lapply(plot_data, function(column) {
      if (is.list(column)) {
        return(column[[i]])
      } else {
        return(column[i])
      }
    }))
    quiet(do.call(plot_func, args))
    grDevices::dev.off()
    if (! file.exists(temp_path)) {
      grDevices::png(temp_path, ...)
      text_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 4, y = 25, size=8, label = "No Plot.") +
        ggplot2::theme_void()
      print(text_plot)
      grDevices::dev.off()
    }
    output <- base64enc::base64encode(temp_path)
    file.remove(temp_path)
    return(paste0('data:image/png;base64,', output))

  }))
  # Dont show selectors with only one option
  if (hide_single_selector) {
    selector <- selector[lapply(selector, length) > 1]
  }
  plot_data[names(selector)] <- lapply(plot_data[names(selector)], as.character)
  plot_data$plot_id <- apply(plot_data[names(selector)], MARGIN = 1, paste0, collapse = '-')

  cat('\n<!--html_preserve-->\n')

  # Make selectors
  selector_class_id <- paste0(imglist_class, '-', names(selector))
  selector_class_id <- gsub(selector_class_id, pattern = '[^a-zA-Z0-9-]+', replacement = '-')
  names(selector_class_id) <- names(selector)
  for (selector_id in names(selector)) {
    cat(paste0('<b>', selector_id, ':  </b>\n'))
    cat(paste0('<select id="', selector_class_id[selector_id], '">\n'))
    cat(paste0('  <option value="', selector[[selector_id]], '">', selector[[selector_id]], '</option>', collapse = '\n'))
    cat(paste0('\n</select>\n'))
  }

  # Add zooming widget
  # TODO: modify so that it is not relying on URLs to stuff on the internet
  if (zoom) {
    cat(paste0('<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/zoomist@2/zoomist.css" />\n'))
    cat(paste0('<script type="module">\n'))
    cat(paste0('  import Zoomist from "https://cdn.jsdelivr.net/npm/zoomist@2/zoomist.js"\n'))
    cat(paste0('  new Zoomist(".zoomist-container-', id_prefix, '", {\n'))
    cat(paste0('    maxScale: 10,\n'))
    cat(paste0('    bounds: true,\n'))
    cat(paste0('    slider: ', tolower(zoom_slider), ',\n'))
    cat(paste0('    zoomer: true\n'))
    cat(paste0('  })\n'))
    cat(paste0('</script>\n'))
    cat(paste0('<div class="zoomist-container-', id_prefix, '">\n'))
    cat(paste0('  <div class="zoomist-wrapper">\n'))
    cat(paste0('    <div class="zoomist-image">\n'))
  }

  # Make image showing the plot
  img_elem_id = paste0('img-', id_prefix)
  cat(paste0('      <img id="', img_elem_id, '" width="100%" src="', plot_data$base64_plot[1], '" />\n'))

  # Add zooming widget
  if (zoom) {
    cat(paste0('    </div>\n'))
    cat(paste0('  </div>\n'))
    cat(paste0('</div>\n'))
  }

  # Add javascript to change with image is shown based on the selector
  function_name <- paste0(gsub(id_prefix, pattern = '[^a-zA-Z0-9]+', replacement = '_'), '_setImgSrc')
  cat(paste0('<script type="text/javascript">\n'))
  cat(paste0('  function ', function_name, '(event) {\n'))
  cat(paste0('    var plots = {', paste0('"', plot_data$plot_id, '": "', plot_data$base64_plot, '"', collapse = ', '), '};\n'))
  cat(paste0('    var selectorIds = [', paste0('"', selector_class_id, '"', collapse = ", "), '];\n'))
  cat(paste0('    var img = document.getElementById("', img_elem_id, '");\n'))
  cat(paste0('    var selectors = selectorIds.map((id) => document.getElementById(id));\n'))
  cat(paste0('    var plot_id = selectors.map((selector) => selector.options[selector.selectedIndex].value).join("-");\n'))
  cat(paste0('    img.src = plots[plot_id];\n'))
  cat(paste0('    return false;\n'))
  cat(paste0('  }\n'))
  cat(paste0('  document.getElementById("', selector_class_id, '").onchange = ', function_name, ';', collapse = '\n'))
  cat(paste0('\n</script>\n'))

  cat('\n<!--/html_preserve-->\n')
}


#' Pick best sendsketch hits for each sample
#'
#' This function processes sendsketch data from pathogen surveillance nextflow pipeline and
#' filters the output so that only the best hits are present
#'
#' @param sketch_data A dataframe containing sketch analysis results.
#' @param sort_columns Character vector; specifies the columns to sort the table by in descending order.
#'   The default is `c("sample_id", "WKID", "ANI", "Complt")`. The `sample_id` is expected to be
#'   unique for each row, and the table will be sorted by this column first, followed by the other
#'   metrics in the order provided.
#' @param top_n Integer; the number of top entries to retain for each unique sample ID after sorting.
#'   The default is `1`, meaning only the top entry is kept.
#'
#' @return A `data.frame`
#'
#' @keywords internal
sendsketch_best_hits <- function(sketch_data, sort_columns = c("WKID", "ANI", "Complt"), top_n = 1) {
  order_data <- c(list(decreasing = TRUE), unname(sketch_data[sort_columns]))
  sketch_data <- sketch_data[do.call(order, order_data), , drop = FALSE]
  split_data <- split(sketch_data, sketch_data[c('sample_id', 'outdir_path')], drop = TRUE)
  final_table <- do.call(rbind, lapply(split_data, function(x) {
    x[seq_len(top_n), , drop = FALSE]
  }))
  rownames(final_table) <- NULL
  return(final_table)
}



#' Combines multiple data frames
#' 
#' combines multiple data frames, handling cases where they have different columns by filling missing values with NA
#'  
#' @param ... One or more data.frames to combine
#'
#' @keywords internal
combine_data_frames <- function(...) {
  # Get all input data frames as a list
  dfs <- list(...)
  
  # Check if all inputs are data frames
  if (!all(sapply(dfs, is.data.frame))) {
    stop("All arguments must be data.frames")
  }
  
  # Get all unique column names across all data frames
  all_cols <- unique(unlist(lapply(dfs, names)))
  
  # Function to align columns of a single data frame
  align_df <- function(df) {
    # Find missing columns in this df
    missing_cols <- setdiff(all_cols, names(df))
    
    # Add missing columns filled with NA
    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        df[[col]] <- NA
      }
    }
    
    # Reorder columns to match all_cols order
    df[, all_cols, drop = FALSE]
  }
  
  # Apply alignment to all data frames and combine them
  combined <- do.call(rbind, lapply(dfs, align_df))
  
  # Reset row names
  rownames(combined) <- NULL
  
  return(combined)
}


#' Combine a list of matrices
#'
#' Combine a list of matrices
#'
#' @param matrices list of matrices
#' @param as_edge_list If `TRUE`, return a table with columns for the row,
#'   column, and value for every unique cell in the combined matrix instead of a
#'   matrix.
#'
#' @keywords internal
combine_matrices <- function(matrices, as_edge_list = FALSE) {
  # Validate input
  if (! is.list(matrices) || ! all(sapply(matrices, is.matrix))) {
    stop("Input must be a list of matrices")
  }
  
  # Convert to edge list
  edges <- do.call(rbind, lapply(matrices, function(x) {
    as.data.frame(as.table(x))
  }))
  colnames(edges) <- c("row", "col", "value")
  all_ids <- unique(c(edges$col, edges$row))
  
  # Check for duplicated row/col with different values
  edges <- unique(edges)
  is_duplicated <- duplicated(edges[, 1:2])
  if (any(is_duplicated)) {
    warning("Shared column names with differing values detected. Removing duplicated values.")
  }
  edges <- edges[! is_duplicated, , drop = FALSE]
  
  if (as_edge_list) {
    return(edges)
  } else {
    # Add missing combinations
    filler <- expand.grid(all_ids, all_ids)
    filler$value <- NA
    colnames(filler) <- c("row", "col", "value")
    edges <- rbind(edges, filler)
    edges <- edges[! duplicated(edges[, 1:2]), , drop = FALSE]
    
    # Combine into a matrix
    edges <- edges[order(match(edges$col, all_ids), match(edges$row, all_ids)), , drop = FALSE]
    combined_matrix <- edges$value
    dim(combined_matrix) <- c(length(all_ids), length(all_ids))
    rownames(combined_matrix) <- all_ids
    colnames(combined_matrix) <- all_ids
    return(combined_matrix)
  }
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
postprocess_table_list <- function(table_list, simplify) {
  for (outdir_path in names(table_list)) {
    if (! 'outdir_path' %in% colnames(table_list[[outdir_path]])) {
      table_list[[outdir_path]]$outdir_path <- outdir_path
    }
  }
  if (simplify) {
    table_list <- do.call(combine_data_frames, table_list)
    col_order <- c(colnames(table_list)[colnames(table_list) != 'outdir_path'], 'outdir_path')
    table_list <- table_list[, col_order, drop = FALSE]
  }
  return(table_list)
}