#' Plot generic phylogeny
#'
#' Plot phylogenies present in the output of a pathogensurveillance run.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output or paths to tree files.
#' @param target The names of one or more output types to search for.
#' @param tree One or more trees to plot. Only used if `path` is not used.
#' @param sample_meta The path to sample metadata or the data itself in the form
#'   of a table. Only used if `path` is not used.
#' @param ref_meta The path to reference metadata or the data itself in the form
#'   of a table. Only used if `path` is not used.
#' @param collapse_by_tax A table with taxonomic data encoded as one column per
#'   rank and a `sample_id` column with the same values present in
#'   `sample_meta`. Used to combine the trees together based on the smallest shared
#'   taxonomic tree. Only used if `path` is not used.
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#'
#' @return  A list of plots, unless `collapse_by_tax` is used, in which case a
#'   single plot is returned.
#'
#' @export
plot_tree <- function(path = NULL, target = NULL, tree = NULL,
                                  sample_meta = NULL, ref_meta = NULL,
                                  collapse_by_tax = NULL, interactive = FALSE) {
  # Find and parse needed data
  if (! is.null(path)) {
    sample_meta <- find_ps_data(path, target = 'sample_metadata', simplify = TRUE)
    ref_meta <- find_ps_data(path, target = 'reference_metadata', simplify = TRUE)
    tree <- find_ps_data(path, target = target, simplify = FALSE)
    tree_path_data <- find_ps_paths(path, target = target, simplify = TRUE)
    
    # Name trees by path metadata (TODO: find more generic approach if possible.)
    if (target == 'multigene_tree') {
      names(tree) <- tree_path_data$method
      if (length(unique(tree_path_data$cluster_id)) > 1) {
        names(tree) <- paste0(names(tree), '_', tree_path_data$cluster_id)
      }
    } else if (target == 'variant_tree') {
      names(tree) <- ref_meta$ref_name[match(tree_path_data$reference_id, ref_meta$ref_id)]
    }
    
  } else {
    # If a single tree is supplied, convert to list
    if (inherits(tree, "phylo")) {
      tree <- list(tree)
    }
  }

  # If no trees are found, return an empty list
  if (length(tree) == 0) {
    return(list())
  }

  # Find which columns are used to provide colors to the trees, if any
  sample_meta_subset <- sample_meta[, 'sample_id', drop = FALSE]
  names(sample_meta_subset) <- 'id'
  sample_meta_subset$name <- sample_meta$name
  if ('description' %in% colnames(sample_meta) && any(sample_meta_subset$name != sample_meta$description)) {
    sample_meta_subset$description <- sample_meta$description
  }  
  ids_in_trees <- unique(unlist(lapply(tree, function(t) t$tip.label)))
  color_by_text <- as.character(sample_meta$color_by[sample_meta$sample_id %in% ids_in_trees])
  color_by_cols <- unique(unlist(strsplit(color_by_text, split = ';')))
  color_by_cols <- color_by_cols[!is.na(color_by_cols)]
  if (length(color_by_cols) > 0) {
    sample_meta_subset <- cbind(sample_meta_subset,  sample_meta[, color_by_cols, drop = FALSE])
  }
  ref_meta_subset <- ref_meta[, c('ref_id', 'ref_name', 'ref_description'), drop = FALSE]
  colnames(ref_meta_subset) <- c('id', 'name', 'description')
  ref_meta_subset <- as.data.frame(do.call(cbind, lapply(colnames(sample_meta_subset), function(col) {
    if (col %in% colnames(ref_meta_subset)) {
      return(ref_meta_subset[[col]])
    } else {
      return(rep(NA, nrow(ref_meta_subset)))
    }
  })))
  colnames(ref_meta_subset) <- colnames(sample_meta_subset)
  ref_meta_subset$`_Sample_type_` <- 'Reference' 
  sample_meta_subset$`_Sample_type_` <- 'Sample'
  metadata <- rbind(sample_meta_subset, ref_meta_subset)
  metadata[] <- lapply(metadata, function(x) {
    x[is.na(x)] <- ''
    x
  })
  
  # Remove the "Root" node label if present 
  tree <- lapply(tree, function(t) {
    t$node.label[t$node.label == "Root"] <- ''
    return(t)
  })
  
  heattree::heat_tree(
    tree = tree,
    metadata = list(metadata)[rep(1, length(tree))],
    aesthetics = list(c('tipLabelColor' = '_Sample_type_', 'tipLabelText' = 'name'))[rep(1, length(tree))],
    manualZoomAndPanEnabled = FALSE
  )
}



#' Plot map widget
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output.
#' @param def_size The default marker size-by variable
#' @param def_color The default marker color-by variable.
#' @param width The width of the widget.
#' @param height The height of the widget.
#' @param elementID Widget element ID.
#'
#' @return  A map plot.
#'
#' @export

plot_map_widget <- function(path = NULL, width = NULL, height = NULL, elementID = NULL, def_size = '', def_color = '') {
  if (!is.null(path)) {
    sample_meta <- find_ps_data(path, target = 'sample_metadata', simplify = TRUE)
    metaMapWidgetR::meta_map_widget(sample_meta, width, height, elementID, def_size, def_color)
  } 
}


#' Plot an interactive MSN of variant data
#'
#' Plot a minimum spanning network for each group of samples aligned to a
#' reference found in `pathogensurveillance` output. The network controls and
#' metadata coloring are provided by `msnwidgetR`.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param interactive Whether to return an interactive HTML widget. Static MSN
#'   plots are not currently supported.
#'
#' @return An interactive `msnwidgetR` widget, or `NULL` when no alignments
#'   contain at least two samples after removing their reference sequence.
#'
#' @examples
#' \dontrun{
#' variant_msn_plot("path/to/pathogensurveillance/output", interactive = TRUE)
#' }
#'
#' @export
variant_msn_plot <- function(path, interactive = TRUE) {
  if (!isTRUE(interactive)) {
    stop("Static MSN plots are not supported; use `interactive = TRUE`.", call. = FALSE)
  }

  alignment_paths <- find_ps_paths(
    path,
    target = "variants_fasta",
    simplify = TRUE,
    must_exist = FALSE
  )
  if (is.null(alignment_paths) || nrow(alignment_paths) == 0L) {
    return(NULL)
  }

  alignments <- find_ps_data(path, target = "variants_fasta")
  sample_data <- find_ps_data(path, target = "sample_metadata", simplify = TRUE)
  reference_data <- find_ps_data(path, target = "reference_metadata", simplify = TRUE)

  datasets <- lapply(seq_len(nrow(alignment_paths)), function(i) {
    alignment_path <- alignment_paths$path[i]
    alignment <- alignments[[alignment_path]]
    if (is.null(alignment) || !length(alignment)) {
      return(NULL)
    }

    reference_id <- as.character(alignment_paths$reference_id[i])
    alignment_ids <- rownames(alignment)
    sample_rows <- !is.na(alignment_ids) & alignment_ids != reference_id
    alignment <- alignment[sample_rows, , drop = FALSE]
    alignment_ids <- alignment_ids[sample_rows]
    if (length(alignment_ids) < 2L) {
      return(NULL)
    }

    alignment <- as.character(alignment)
    sequences <- apply(alignment, 1L, paste0, collapse = "")

    metadata <- data.frame(id = alignment_ids, stringsAsFactors = FALSE)
    if (!is.null(sample_data) && "sample_id" %in% names(sample_data)) {
      sample_rows <- match(alignment_ids, as.character(sample_data$sample_id))
      metadata <- as.data.frame(sample_data[sample_rows, , drop = FALSE])
      metadata$sample_id <- NULL
      metadata$id <- alignment_ids
      metadata <- metadata[, c("id", setdiff(names(metadata), "id")), drop = FALSE]
    }

    reference_name <- if (
      !is.null(reference_data) &&
        all(c("ref_id", "ref_name") %in% names(reference_data))
    ) {
      reference_data$ref_name[match(reference_id, reference_data$ref_id)]
    } else {
      NA_character_
    }
    if (is.na(reference_name) || !nzchar(reference_name)) {
      reference_name <- reference_id
    }

    list(
      name = paste0(reference_name, " (", reference_id, ")"),
      fasta = paste(paste0(">", alignment_ids, "\n", sequences), collapse = "\n"),
      metadata = metadata
    )
  })
  datasets <- Filter(Negate(is.null), datasets)
  if (!length(datasets)) {
    return(NULL)
  }

  msnwidgetR::msnwidget(datasets = datasets)
}


#' Plot ANI matrix
#'
#' Plot ANI matrix with dendrogram from data present in the output of a
#' pathogensurveillance run.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output or paths to tree files.
#' @param combine If `TRUE`, combine data from all ANI matrices found into a
#'   single plot.
#' @param interactive Whether to return an HTML-based interactive widget
#'   (default: FALSE)
#' @param max_size The maximum number of rows/columns to plot. To subset the
#'   data to this amount when needed, references most similar to the samples
#'   will be used. If there are more samples than this, then only samples will
#'   be plotted. This does not subset samples. See 'prefer_contextual` as well
#' @param prefer_contextual How subsets of rows/columns are chosen when
#'   `max_size` takes effect. If `TRUE` and contextual references have been
#'   chosen, give preference to including them.
#' @param height The height in pixels. If not `interactive`, this is divided by
#'   `dpi` to convert it to inches.
#' @param width The width in pixels. If not `interactive`, this is divided by
#'   `dpi` to convert it to inches.
#' @param dpi How pixels are converted to inches
#' @param font_size Size of text used for labels
#' @param max_label_length Labels longer than this length will be shortened.
#'
#' @return An interactive HTML widget when `interactive = TRUE`, otherwise a
#'   static heatmap.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' estimated_ani_heatmap(path)
#' estimated_ani_heatmap(path, interactive = TRUE)
#'
#' @export
estimated_ani_heatmap <- function(path, combine = FALSE, interactive = FALSE, max_size = 100,
                                  prefer_contextual = TRUE, height = NULL, width = NULL,
                                  dpi = 100, font_size = 8, max_label_length = 30) {
  if (combine) {
    stop('The `combine` option is not yet supported.')
  }
  
  # Find and parse data
  input_matrix <- find_ps_data(path, target = 'ani_matrix_csv', simplify = TRUE)
  sample_meta <- find_ps_data(path, target = 'sample_metadata', simplify = TRUE)
  ref_meta <- find_ps_data(path, target = 'reference_metadata', simplify = TRUE)
  sample_ids <- colnames(input_matrix)[colnames(input_matrix) %in% sample_meta$sample_id]
  ref_ids <- colnames(input_matrix)[colnames(input_matrix) %in% ref_meta$ref_id]
  
  if (!is.null(max_size)) {
    if (nrow(input_matrix) > max_size) {
      refs_to_keep <- min(max_size - length(sample_ids), length(ref_ids))
      if (refs_to_keep > 0) {
        # Calculate mean similarity for each reference to all samples
        ref_similarities <- vapply(ref_ids, FUN.VALUE = numeric(1), function(ref_id) {
          similarities <- input_matrix[sample_ids, ref_id, drop = TRUE]
          mean(similarities, na.rm = TRUE)
        })
        
        # Check if contextual references can/should be preferred
        ref_path_data <- find_ps_data(path, target = 'contextual_refs', simplify = TRUE)
        if (prefer_contextual && !is.null(ref_path_data)) {
          is_contextual <- ref_ids %in% ref_path_data$reference_id
          ordered_refs <- ref_ids[order(is_contextual, ref_similarities, decreasing = TRUE, na.last = TRUE)]
        } else {
          ordered_refs <- ref_ids[order(ref_similarities, decreasing = TRUE, na.last = TRUE)]
        }
        
        # Subset to best refs
        ids_used <- c(sample_ids, ordered_refs[1:refs_to_keep])
        input_matrix <- input_matrix[ids_used, ids_used]
      } else { # If we have too many samples, just use samples
        input_matrix <- input_matrix[sample_ids, sample_ids]
      }
    }
  }
  
  make_heatmap(input_matrix = input_matrix, sample_data = sample_meta, ref_data = ref_meta,
               interactive = interactive, height = height, width = width, dpi = dpi,
               font_size = font_size, max_label_length = max_label_length,
               metric = "ani")
}


#' Plot POCP matrix
#'
#' Plot POCP matrix with dendrogram from data present in the output of a
#' pathogensurveillance run.
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output or paths to tree files.
#' @param combine If `TRUE`, combine data from all POCP matrices found into a
#'   single plot.
#' @param interactive Whether to return an HTML-based interactive widget
#'   (default: FALSE)
#' @param height The height in pixels. If not `interactive`, this is divided by
#'   `dpi` to convert it to inches.
#' @param width The width in pixels. If not `interactive`, this is divided by
#'   `dpi` to convert it to inches.
#' @param dpi How pixels are converted to inches
#' @param font_size Size of text used for labels
#' @param max_label_length Labels longer than this length will be shortened.
#'
#' @return An interactive HTML widget when `interactive = TRUE`, otherwise a
#'   static heatmap.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' pocp_heatmap(path)
#' pocp_heatmap(path, interactive = TRUE)
#'
#' @export
pocp_heatmap <- function(path, combine = FALSE, interactive = FALSE,
                         height = NULL, width = NULL, dpi = 100, font_size = 8,
                         max_label_length = 30) {
  if (combine) {
    stop('The `combine` option is not yet supported.')
  }

  # Find and parse data
  input_matrix <- find_ps_data(path, target = 'pocp', simplify = TRUE)
  sample_meta <- find_ps_data(path, target = 'sample_metadata', simplify = TRUE)
  ref_meta <- find_ps_data(path, target = 'reference_metadata', simplify = TRUE)
  
  make_heatmap(input_matrix = input_matrix, sample_data = sample_meta, ref_data = ref_meta,
               interactive = interactive, height = height, width = width, dpi = dpi,
               font_size = font_size, max_label_length = max_label_length,
               metric = "pocp")
}



#' Make an heatmap with dendrogram.
#'
#' The heatmap is based on a distance matrix
#'
#' @param input_matrix A square similarity matrix
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param ref_data A tibble/data.frame with information on references used in analysis
#' @param interactive Whether to return an HTML-based interactive widget (default: FALSE)
#' @param height The height in pixels. If not `interactive`, this is divided by `dpi` to convert it to inches.
#' @param width The width in pixels. If not `interactive`, this is divided by `dpi` to convert it to inches.
#' @param dpi How pixels are converted to inches
#' @param font_size Size of text used for labels
#' @param max_label_length Labels longer than this length will be shortened for display.
#' @param metric Similarity metric used to label the interactive widget.
#'
#' @return An interactive HTML widget when `interactive = TRUE`, otherwise a
#'   static heatmap.
#'
#' @keywords internal
make_heatmap <- function(input_matrix, ref_data, sample_data, interactive = FALSE,
                         height = NULL, width = NULL, dpi = 100, font_size = 10,
                         max_label_length = 30, metric = c("ani", "pocp")) {
  metric <- match.arg(metric)
  input_matrix <- as.matrix(input_matrix)
  storage.mode(input_matrix) <- "numeric"
  observed_values <- input_matrix[is.finite(input_matrix)]
  if (metric == "ani" && length(observed_values) && max(observed_values) <= 1) {
    input_matrix <- input_matrix * 100
  }
  if (nrow(input_matrix) != ncol(input_matrix)) {
    stop("The similarity matrix must be square.", call. = FALSE)
  }
  if (is.null(rownames(input_matrix)) || is.null(colnames(input_matrix))) {
    stop("The similarity matrix must have row and column names.", call. = FALSE)
  }
  if (!setequal(rownames(input_matrix), colnames(input_matrix))) {
    stop("The matrix row and column names must contain the same identifiers.", call. = FALSE)
  }
  if (!identical(rownames(input_matrix), colnames(input_matrix))) {
    input_matrix <- input_matrix[rownames(input_matrix), rownames(input_matrix), drop = FALSE]
  }

  name_key <- c(
    if (!is.null(ref_data) && nrow(ref_data) && all(c("ref_id", "ref_name") %in% names(ref_data))) {
      stats::setNames(as.character(ref_data$ref_name), as.character(ref_data$ref_id))
    },
    if (!is.null(sample_data) && nrow(sample_data) && all(c("sample_id", "name") %in% names(sample_data))) {
      stats::setNames(as.character(sample_data$name), as.character(sample_data$sample_id))
    }
  )
  ids <- rownames(input_matrix)
  full_labels <- if (length(name_key)) {
    unname(name_key[ids])
  } else {
    rep(NA_character_, length(ids))
  }
  missing_names <- is.na(full_labels) | full_labels == ""
  full_labels[missing_names] <- ids[missing_names]
  duplicated_names <- duplicated(full_labels) | duplicated(full_labels, fromLast = TRUE)
  full_labels[duplicated_names] <- paste0(full_labels[duplicated_names], " (", ids[duplicated_names], ")")
  label_limit <- max(1L, as.integer(max_label_length))
  display_labels <- ifelse(
    nchar(full_labels) > label_limit,
    paste0(substr(full_labels, 1L, label_limit), "\u2026"),
    full_labels
  )
  display_labels <- make.unique(display_labels, sep = " ")
  sample_ids <- if (!is.null(sample_data) && "sample_id" %in% names(sample_data)) {
    intersect(ids, as.character(sample_data$sample_id))
  } else {
    character()
  }
  if (interactive) {
    metric_label <- if (metric == "ani") "Estimated ANI similarity (%)" else "Percentage of conserved proteins (%)"
    matrix_payload <- lapply(seq_len(nrow(input_matrix)), function(row) {
      lapply(seq_len(ncol(input_matrix)), function(column) {
        value <- input_matrix[row, column]
        if (length(value) == 0L || is.na(value) || !is.finite(value)) NULL else unname(value)
      })
    })
    return(htmlwidgets::createWidget(
      name = "psHeatmap",
      x = list(
        labels = ids,
        display_labels = unname(display_labels),
        full_labels = unname(full_labels),
        sample_ids = unname(sample_ids),
        matrix = matrix_payload,
        metric = metric,
        metric_label = metric_label,
        font_size = min(14, max(4, as.numeric(font_size)))
      ),
      width = width,
      height = height,
      package = "PathoSurveilR",
      sizingPolicy = htmlwidgets::sizingPolicy(
        defaultWidth = if (is.null(width)) 950 else width,
        defaultHeight = if (is.null(height)) 850 else height,
        viewer.fill = TRUE,
        browser.fill = TRUE,
        padding = 0
      )
    ))
  }

  width_inches <- if (is.null(width)) NA_real_ else width / dpi
  height_inches <- if (is.null(height)) NA_real_ else height / dpi
  na_to_zero <- input_matrix
  na_to_zero[!is.finite(na_to_zero)] <- NA_real_
  na_to_zero[is.na(na_to_zero)] <- 0
  dimnames(input_matrix) <- list(display_labels, display_labels)
  stats::hclust(stats::dist(na_to_zero), method = "complete") -> clustered
  pheatmap::pheatmap(input_matrix,
                     cluster_rows = clustered,
                     cluster_cols = clustered,
                     angle_col = 45,
                     show_rownames = TRUE, labels_row = display_labels,
                     labels_col = display_labels,
                     fontsize_row = font_size, fontsize_col = font_size,
                     width = width_inches, height = height_inches)
}


#' Make sunburst plot of sendsketch taxonomy
#'
#' Converts classifications of top hits in sendsketch output into an interactive
#' sunburst plot.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output or a table in the format of the
#'   [sendsketch_parsed()] output.
#' @param interactive Whether or not to produce an interactive
#'   HTML/javascript-based plot or a static one.
#' @param ... Passed to `sendsketch_best_hits`
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' sendsketch_taxonomy_plot(path)
#' sendsketch_taxonomy_plot(path, interactive = TRUE)
#' \dontshow{
#' ## R CMD check: make sure any open connections are closed afterward
#' closeAllConnections()
#' }
#' @export
sendsketch_taxonomy_plot <- function(path, interactive = FALSE, ...) {

  # Parse the input if it is a file/folder path
  if (is.data.frame(path)) {
    sketch_data <- path
  } else {
    sketch_data <- sendsketch_parsed(path, only_best = TRUE)
  }

  # Sort and filter data
  top_hits <- sendsketch_best_hits(sketch_data, ...)

  # Make table with values duplicated for each taxon in all classifications
  split_taxa <- strsplit(top_hits$taxonomy, split = ';')
  tax_data <- do.call(rbind, lapply(seq_len(nrow(top_hits)), function(i) {
    n_taxa <- length(split_taxa[[i]])
    out <- top_hits[rep(i, n_taxa), , drop = FALSE]
    out$taxon <- vapply(seq_len(n_taxa), FUN.VALUE = character(1), function(n) {
      paste0(split_taxa[[i]][1:n], collapse = ';')
    })
    out$supertaxon <- c(NA_character_, out$taxon[1:(n_taxa - 1)])
    out$tip <- vapply(seq_len(n_taxa), FUN.VALUE = character(1), function(n) {
      split_taxa[[i]][n]
    })
    rownames(out) <- NULL
    return(out)
  }))

  # Convert to edge list
  plot_data <- unique(tax_data[, c('supertaxon', 'taxon', 'tip'), drop = FALSE])
  plot_data$count <- vapply(plot_data$taxon, FUN.VALUE = numeric(1), function(t) {
    sum(t == tax_data$taxon)
  })

  # Make unique IDs for each taxon
  id_key <- as.character(seq_len(nrow(plot_data)))
  names(id_key) <- plot_data$taxon
  plot_data$from <- id_key[plot_data$supertaxon]
  plot_data$to <- id_key[plot_data$taxon]

  # Rename taxon to just the tip taxon name and include rank if needed to make unique
  plot_data$rank <- sub(plot_data$tip, pattern = '^([a-z]+?):(.+)$', replacement = '\\1')
  plot_data$name <- sub(plot_data$tip, pattern = '^([a-z]+?):(.+)$', replacement = '\\2')
  plot_data$is_unique <- vapply(plot_data$name, FUN.VALUE = logical(1), function(x) {
    sum(x == plot_data$name) <= 1
  })
  plot_data$name <- ifelse(plot_data$is_unique, plot_data$name, plot_data$tip)

  output <- plotly::plot_ly(
    type = 'sunburst',
    ids = plot_data$to,
    labels = plot_data$name,
    parents = plot_data$from,
    values = plot_data$count,
    domain = list(column = 0),
    branchvalues = 'total'
  )

  if (! interactive) {
    temp_file_html <- tempfile(fileext = ".html")
    temp_file_png <- tempfile(fileext = ".png")
    htmlwidgets::saveWidget(widget = plotly::config(output, displayModeBar = FALSE), file = temp_file_html)
    output <- webshot2::webshot(url = temp_file_html, file = temp_file_png,
                                delay = 1, vheight = 750, vwidth = 750, zoom = 2)
  }

  return(output)
}
