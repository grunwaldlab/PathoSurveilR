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
    tree <- find_ps_data(path, target = target)
    collapse_by_tax <- sendsketch_parsed(path, only_best = TRUE, only_shared = TRUE, update_taxonomy = TRUE, simplify = TRUE)
    rank_cols <- colnames(collapse_by_tax)[39:ncol(collapse_by_tax)]
    collapse_by_tax <- collapse_by_tax[, c('sample_id', rank_cols), drop = FALSE]
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
  ids_in_trees <- unique(unlist(lapply(tree, function(t) t$tip.label)))
  color_by_text <- as.character(sample_meta$color_by[sample_meta$sample_id %in% ids_in_trees])
  color_by_cols <- unique(unlist(strsplit(color_by_text, split = ';')))
  color_by_cols <- color_by_cols[!is.na(color_by_cols)]
  color_by_col_names <- c(color_by_cols, 'Default')
  color_by_cols <- c(as.list(color_by_cols), list(NULL))  # NULL ensures that the default color scheme is also used

  # Plot one tree for each color_by column
  tree_plots <- lapply(color_by_cols, function(color_by) {
    plot_phylogeny(
      tree,
      sample_meta,
      ref_meta,
      color_by,
      collapse_by_tax, 
      interactive = interactive
    )
  })
  names(tree_plots) <- color_by_col_names
  return(tree_plots)
}


#' Plot a phylogeny
#'
#' @param trees One or more [ape::phylo()] objects. Inputting multiple trees
#'   only makes sense when `collapse_by_tax` is used.
#' @param sample_meta A table containing metadata for samples in `tree_paths`.
#' @param ref_meta A table containing metadata for references in `tree_paths`.
#' @param collapse_by_tax A [tibble::tibble()] with the columns `sample_id` and
#'   those named by taxonomic ranks. These are ids present in `sample_meta`.
#'   These are used to provide a taxonomic tree that functions as a backbone to
#'   combine multiple trees into a single one. (Default: return a list of plots
#'   instead of a single plot)
#' @param color_by The name of the column in the metadata used to color samples.
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#' @param tip_label_size The size of tree tip labels. 1-5 might be reasonable values.
#' @param node_label_size The size of tree node (bootstrap) labels. 1-5 might be reasonable values.
#'
#' @import ggtree
#'
#' @keywords internal
plot_phylogeny <- function(trees, sample_meta, ref_meta, color_by = NULL, collapse_by_tax = NULL,
                           interactive = FALSE, tip_label_size = NULL, node_label_size = NULL, ...) {
  if (interactive) {
    stop('Interactive tree plotting is not yet supported.')
  }

  # If a single tree is provided, convert it to a list of a single tree to simplify code below
  if (inherits(trees, "phylo")) {
    trees <- list(trees)
  }
  
  # Combine trees by connecting to the taxonomy-derived tree
  if (is.null(collapse_by_tax) || length(trees) == 1) {
    if (length(trees) == 1) {
      combined_tree <- trees[[1]]
    } else {
      stop(call. = FALSE, 'If "collapse_by_tax" is not provided, then only a single tree can be plotted at a time. Multiple trees were given.')
    }
  } else {

    # If taxonomy does not have at least one conserved rank, then add one
    colnames(collapse_by_tax) <- gsub(pattern = ' ', replacement = '_', colnames(collapse_by_tax))
    ranks <- colnames(collapse_by_tax)[colnames(collapse_by_tax) != 'sample_id']
    n_unique_taxa <- unlist(lapply(collapse_by_tax[ranks], function(x) length(unique(x))))
    if (! any(n_unique_taxa == 1)) {
      collapse_by_tax <- cbind(sample_id = collapse_by_tax$sample_id, root = 'Life', collapse_by_tax[ranks])
      ranks <- c('root', ranks)
    }

    # Subset taxonomy data to just the part used by each tree
    tree_tax_data <- do.call(rbind, lapply(seq_along(trees), function(i) {
      tree_ids <- trees[[i]]$tip.label
      tree_samp_ids <- tree_ids[tree_ids %in% collapse_by_tax$sample_id]
      out <- cbind(tree = i, collapse_by_tax[match(tree_samp_ids, collapse_by_tax$sample_id), ranks])
      is_variable_col <- vapply(out[ranks], FUN.VALUE = logical(1), function(y) length(unique(y)) > 1)
      variable_cols <- ranks[is_variable_col]
      last_uniform_taxon <- ranks[max(which(! is_variable_col))]
      out[, variable_cols] <- unique(out[, last_uniform_taxon])
      unique(out)
    }))
    tree_tax_data[ranks] <- lapply(tree_tax_data[ranks], as.factor)

    # Combine trees
    base_tree <- ape::as.phylo(stats::as.formula(paste0('~', paste0(ranks, collapse = '/'))), data = tree_tax_data)
    base_tree$tip.label <- make.unique(base_tree$tip.label)
    base_tree$node.label <- make.unique(base_tree$node.label)
    max_edge_len <- max(unlist(lapply(trees, function(x) ape::node.depth.edgelength(x)[1:ape::Ntip(x)])))
    base_tree$edge.length <- rep(max_edge_len * 0.1, nrow(base_tree$edge))
    combined_tree <- base_tree
    tip_rank <- as.character(tree_tax_data[, ranks[length(ranks)]])
    index_key <- match(base_tree$tip.label, make.unique(tip_rank))
    for (index in rev(seq_along(trees))) {
      combined_tree <- ape::bind.tree(combined_tree, trees[[index_key[index]]], where = index)
    }
  }

  # Prepare labels
  label_key <- c(
    stats::setNames(sample_meta$description, sample_meta$sample_id),
    stats::setNames(ref_meta$ref_description, ref_meta$ref_id)
  )
  tip_labels <- label_key[combined_tree$tip.label]
  tip_labels <- stats::setNames(make.unique(tip_labels, sep = ' '), names(tip_labels))

  # Prepare colors
  if (is.null(color_by)) {
    color_by <- '_sequence_type_'
    sample_meta[['_sequence_type_']] <- 'Sample'
    ref_meta[['_sequence_type_']] <- 'Reference'
  }
  ids_in_trees <- combined_tree$tip.label
  color_key <- character(0)
  if (color_by %in% colnames(sample_meta)) {
    color_key <- c(color_key, stats::setNames(sample_meta[[color_by]], sample_meta$sample_id))
  }
  if (color_by %in% colnames(ref_meta)) {
    color_key <- c(color_key, stats::setNames(ref_meta[[color_by]], ref_meta$ref_id))
  }
  tip_colors <- color_key[combined_tree$tip.label]

  # Make data associated with tree tips
  tip_data <- tibble::tibble(
    newick_label = combined_tree$tip.label,
    tip_color = unname(tip_colors),
    tip_label = tip_labels[combined_tree$tip.label]
  )

  # Get tip labels from input trees to find their MRCA nodes in combined tree
  input_tree_tip_groups <- list()
  if (length(trees) > 1) {
    input_tree_tip_groups <- lapply(trees, function(tree) {
      tree$tip.label
    })
  }

  if (is.null(combined_tree$node.label)) {
    all_labels = c(combined_tree$tip.label, rep(NA, combined_tree$Nnode))
  } else {
    all_labels = c(combined_tree$tip.label, combined_tree$node.label)
  }

  # Use groupOTU with tip labels to identify branches that belong to input trees
  if (length(input_tree_tip_groups) > 0) {
    # Create groups based on tip labels from each input tree
    input_tree_groups <- ggtree::groupOTU(combined_tree, input_tree_tip_groups)
    group_info <- attr(input_tree_groups, "group")
    
    # Branches with input tree tips should be black/solid, base tree should be grey/dashed
    is_input_tree_branch <- group_info > 0
    
    node_data <- tibble::tibble(
      newick_label = all_labels,
      node_label = ifelse(is.na(all_labels), "", all_labels),
      branch_color = ifelse(is_input_tree_branch, "black", "grey"),
      branch_type = ifelse(is_input_tree_branch, "solid", "dashed")
    )
  } else {
    # Single tree case - all branches black and solid
    node_data <- tibble::tibble(
      newick_label = all_labels,
      node_label = ifelse(is.na(all_labels), "", all_labels),
      branch_color = "black",
      branch_type = "solid"
    )
  }

  legend_title <- tools::toTitleCase(trimws(gsub(color_by, pattern = '_', replacement = ' ')))

  logistic_scaling_func <- function(shared_count, floor = 0.1, ceiling = 2.5, midpoint = 0, steepness = 0.02) {
    logistic_value <- ceiling / (1 + exp(1)^(-steepness * (shared_count - midpoint)))
    ((logistic_value - 0.5 + floor) / (0.5 - floor)) - ceiling + 1
  }

  if (is.null(node_label_size)) {
    node_label_size = 3 - logistic_scaling_func(nrow(tip_data), steepness = 0.005, floor = 0, ceiling = 2.5)
  }
  if (is.null(tip_label_size)) {
    tip_label_size = 4 - logistic_scaling_func(nrow(tip_data), steepness = 0.005, floor = 0, ceiling = 3.5)
  }

  plotted_tree <- ggtree(combined_tree, ggplot2::aes_string(color = 'branch_color', linetype = 'branch_type')) %<+% node_data +
    geom_nodelab(ggplot2::aes_string(label = 'node_label'), hjust = 1.3, nudge_y = 0.3, size = node_label_size) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_linetype_identity()
  plotted_tree <- plotted_tree %<+% tip_data +
    scale_x_continuous(expand = ggplot2::expansion(add = c(0.1, 0.2 + max(nchar(tip_data$tip_label)) * 0.12 / sqrt(nrow(tip_data))))) +
    ggnewscale::new_scale_color() +
    geom_tiplab(ggplot2::aes_string(label = 'tip_label', color = 'tip_color'), size = tip_label_size) +
    geom_tippoint(ggplot2::aes_string(color = 'tip_color'), alpha = 0) + # Invisible tips just there to make override.aes below change the legend color shapes
    ggplot2::scale_color_viridis_d(end = 0.8, na.value = "black") +
    ggplot2::guides(color = guide_legend(title = legend_title, override.aes = list(label = "", size = 3, alpha = 1, shape = 15))) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(1,1,1,1, "cm")
    )

  return(plotted_tree)
}


#' Plot MSN of variant data
#'
#' Plot a minimum spanning network for each group of samples aligned to a
#' reference found in `pathogensurveillance` output.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output.
#' @param combine If `TRUE` combine multiple MSNs into a single figure and
#'   return a single figure. If `FALSE`, return a list of figures.
#'
#' @examples
#' path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
#' variant_msn_plot(path)
#'
#' @export
variant_msn_plot <- function(path, combine = TRUE) {
  align_data <- find_ps_paths(path, target = 'variants_fasta')
  alignments <- find_ps_data(path, target = 'variants_fasta')
  sample_data <- find_ps_data(path, target = 'sample_metadata', simplify = TRUE)
  ref_data <- find_ps_data(path, target = 'reference_metadata', simplify = TRUE)

  # Find which columns are used to provide colors to the trees, if any
  ids_used <- unique(unlist(lapply(alignments, function(a) {
    if (is.null(a)) {
      return(character(0))
    } else {
      return(rownames(a))
    }
  })))
  color_by_text <- as.character(sample_data$color_by[sample_data$sample_id %in% ids_used])
  color_by_cols <- unique(unlist(strsplit(color_by_text, split = ';')))
  color_by_cols <- color_by_cols[! is.na(color_by_cols)]
  color_by_col_names <- c(color_by_cols, 'Default')
  color_by_cols <- c(as.list(color_by_cols), list(NULL))  # NULL ensures that the default color scheme is also used

  # Plot MSNs
  graphics::plot.new()
  output <- lapply(seq_len(nrow(align_data)), function(i) {
    align_without_ref <- alignments[[i]][rownames(alignments[[i]]) != align_data$reference_id[i], ]
    if (is.null(align_without_ref)) {
      return(NULL)
    }
    plot_data <- make_MSN(align_without_ref, sample_data, user_seed = 1, snp_diff_prop = 0.1, population = NULL)
    if (is.null(plot_data)) {
      return(NULL)
    }
    out <- grDevices::recordPlot()
    graphics::plot.new()
    return(out)
  })

  return(output)
}


#' Make Minimum spanning network
#'
#' @param snp_fasta_alignment A DNA alignment in fasta format.
#' @param sample_data A data frame containing information about samples.
#' @param population A character string specifying the column name in sample_data to be used for stratification.
#' @param interactive A logical value indicating whether Whether or not to produce an interactive HTML/javascript-based figures or tables or a static ones.
#' @param snp_threshold An integer specifying the number of SNPs to be used as the threshold for filtering. User can specify whole numbers or alternatively, a relative proportion using 'snp_diff_prop', but not both. Default is NULL.
#' @param snp_diff_prop A numeric value specifying the proportion of SNPs to be used as the threshold for filtering. If user prefers to specify whole numbers, use 'snp_threshold' instead. Default is NULL.
#' @param use_cutoff_predictor A logical value indicating whether to use cutoff predictor for determining SNP threshold.
#' @param show_MLG_table A logical value indicating whether to display a table of multi-locus genotypes.
#' @param user_seed An optional integer specifying the seed for reproducibility.
#' @param ... Additional arguments to be passed to internal functions.
#' @return Minimum spanning network
#'
#' @keywords internal
make_MSN <- function(snp_fasta_alignment, sample_data, population = NULL, interactive = FALSE, snp_threshold = NULL, snp_diff_prop = NULL, use_cutoff_predictor = FALSE, show_MLG_table = FALSE, user_seed = NULL, ...) {

  set.seed(user_seed)
  snp_aln.gi <- adegenet::DNAbin2genind(snp_fasta_alignment)
  if (is.null(snp_aln.gi)) {
    return(NULL)
  }
  snp_genclone <- poppr::as.genclone(snp_aln.gi)

  if (use_cutoff_predictor) {
    snpdist_stats <- poppr::filter_stats(snp_genclone)
    average_thresh <- poppr::cutoff_predictor(snpdist_stats$average$THRESHOLDS)
    poppr::mlg.filter(snp_genclone, distance = poppr::bitwise.dist, percent = TRUE) <- average_thresh
  } else if (!is.null(snp_threshold)) {
    poppr::mlg.filter(snp_genclone, distance = poppr::bitwise.dist, percent = FALSE) <- snp_threshold
  } else if (!is.null(snp_diff_prop)) {
    poppr::mlg.filter(snp_genclone, distance = poppr::bitwise.dist, percent = TRUE) <- snp_diff_prop
  }

  if (is.null(population)) {
    sample_data$no_factor_provided <- 'Sample (No factor provided)'
    population <- 'no_factor_provided'
  }

  sample_data <- sample_data[sample_data$sample_id %in% rownames(snp_fasta_alignment), , drop = FALSE]
  user_factor <- sample_data[[population]]
  node_color <- as.factor(ifelse(is.na(user_factor) | user_factor == "", "Unknown", user_factor))
  myColors <- grDevices::rainbow(length(unique(node_color)))
  names(myColors) <- levels(node_color)
  num_columns <- ncol(sample_data)
  adegenet::strata(snp_genclone) <- cbind(sample_data[, c(1:num_columns)], color_node_by = node_color)
  adegenet::setPop(snp_genclone) <- ~color_node_by

  if (!is.null(snp_threshold)) {
    ms.loc <- poppr::poppr.msn(snp_genclone,
                               distmat = poppr::bitwise.dist(snp_genclone, percent = FALSE),
                               include.ties = TRUE,
                               showplot = FALSE)
  }
  else {
    ms.loc <- poppr::poppr.msn(snp_genclone,
                               distmat = poppr::bitwise.dist(snp_genclone, percent = TRUE),
                               include.ties = TRUE,
                               showplot = FALSE)
  }

  the_edges <- igraph::E(ms.loc$graph)$weight
  edges <- as.list(the_edges)

  if (length(igraph::V(ms.loc$graph)) > 1) {
    output <- poppr::plot_poppr_msn(
      snp_genclone,
      poppr_msn = ms.loc,
      palette = myColors,
      mlg = FALSE,
      quantiles = FALSE,
      wscale = FALSE,
      inds = "None",
      ...
    )
  } else {
    output <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 25, size=8, label = "All samples are in the same multilocus genotype.") +
      ggplot2::theme_void()
  }


  # if (show_MLG_table) {
  #   idlist <- mlg.id(snp_genclone)
  #   mlglist <- data.frame("MLG", "strain")
  #   colnames(mlglist) <- c("V1", "V2")
  #
  #   for (name in names(idlist)) {
  #     newframe <- as.data.frame(cbind(paste0("MLG", "_", name), idlist[[name]]))
  #     mlglist <- rbind(mlglist, newframe)
  #   }
  #
  #   colnames(mlglist) <- c("Multi-locus genotype", "Strain")
  #   mlglist <- mlglist[mlglist$Strain != "strain",]
  #
  #   if (interactive) {
  #     DT::datatable(mlglist, class = "display nowrap", ...) %>%
  #       formatStyle(colnames(mlglist), "white-space" = "nowrap")
  #
  #   } else {
  #     print(mlglist)
  #   }
  # }

  return(output)
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
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
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
#' @return  A list of plots, unless `combine` is used, in which case a single
#'   plot is returned.
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
               font_size = font_size, max_label_length = max_label_length)
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
#' @param interactive Whether to use an HTML-based interactive format or not
#'   (default: TRUE)
#' @param height The height in pixels. If not `interactive`, this is divided by
#'   `dpi` to convert it to inches.
#' @param width The width in pixels. If not `interactive`, this is divided by
#'   `dpi` to convert it to inches.
#' @param dpi How pixels are converted to inches
#' @param font_size Size of text used for labels
#' @param max_label_length Labels longer than this length will be shortened.
#'
#' @return  A list of plots, unless `combine` is used, in which case a single
#'   plot is returned.
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
  input_matrix <- find_ps_data(path, target = 'ani_matrix_csv', simplify = TRUE)
  sample_meta <- find_ps_data(path, target = 'sample_metadata', simplify = TRUE)
  ref_meta <- find_ps_data(path, target = 'reference_metadata', simplify = TRUE)
  
  make_heatmap(input_matrix = input_matrix, sample_data = sample_meta, ref_data = ref_meta,
               interactive = interactive, height = height, width = width, dpi = dpi,
               font_size = font_size, max_label_length = max_label_length)
}



#' Make an heatmap with dendrogram.
#'
#' The heatmap is based on a distance matrix
#'
#' @param input_matrix Approximate ANI matrix output from Sourmash analysis
#' @param sample_data A tibble/data.frame with the sample metadata
#' @param ref_data A tibble/data.frame with information on references used in analysis
#' @param interactive Whether to use an HTML-based interactive format or not (default: TRUE)
#' @param height The height in pixels. If not `interactive`, this is divided by `dpi` to convert it to inches.
#' @param width The width in pixels. If not `interactive`, this is divided by `dpi` to convert it to inches.
#' @param dpi How pixels are converted to inches
#' @param font_size Size of text used for labels
#' @param max_label_length Labels longer than this length will be shortened.
#'
#' @return A heatmap with dendrogram
#'
#' @keywords internal
make_heatmap <- function(input_matrix, ref_data, sample_data, interactive = FALSE,
                         height = NULL, width = NULL, dpi = 100, font_size = 10,
                         max_label_length = 30) {
  # Rename rows/columns for plotting
  name_key <- c(
    stats::setNames(ref_data$ref_name, ref_data$ref_id),
    stats::setNames(sample_data$name, sample_data$sample_id)
  )
  name_key <- ifelse(nchar(name_key) > max_label_length,
                     paste0(substr(name_key, start = 1, stop = max_label_length), '\u2026'),
                     name_key)
  name_key <- stats::setNames(make.unique(name_key, sep = ' '), names(name_key))
  colnames(input_matrix) <- name_key[colnames(input_matrix)]
  rownames(input_matrix) <- name_key[rownames(input_matrix)]
  na_to_zero <- input_matrix
  na_to_zero[is.na(na_to_zero)] <- 0
  clustered <- stats::hclust(stats::dist(na_to_zero), method = "complete")
  if (interactive) {
    dist_func <- function(x) {
      return(clustered)
    }
    output <- heatmaply::heatmaply(input_matrix,
                                   fontsize_row = font_size, fontsize_col = font_size,
                                   width = width, height = height,
                                   hclustfun = function(x) {clustered},
                                   distfun  = function(x) {stats::dist(na_to_zero)},
                                   grid_color = '#EEEEEE')
  } else {
    if (is.null(width)) {
      width <- NA
    } else {
      width <- width / dpi
    }
    if (is.null(height)) {
      height <- NA
    } else {
      height <- height / dpi
    }
    output <- pheatmap::pheatmap(input_matrix,
                                 cluster_rows = clustered,
                                 cluster_cols = clustered,
                                 angle_col = 45,
                                 show_rownames = TRUE, labels_row = colnames(input_matrix),
                                 fontsize_row = font_size, fontsize_col = font_size,
                                 width = width, height = height)
  }
  return(output)
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

#' Make map widget
#'
#' Converts classifications of top hits in sendsketch output into an interactive
#' sunburst plot.
#'
#' @param path The path to one or more folders that contain
#'   pathogensurveillance output or a sample metadata table with coordinates.
#'
#' @export
sample_distribution_map <- function(path) {
  
  # Parse the input if it is a file/folder path
  # if (is.data.frame(path)) {
  #   metadata <- path
  # } else {
  #   metadata <- find_ps_data(path, target = 'sample_metadata') 
  # }
  
  # TODO: temp test data, remove later
  metadata <- data.frame(
    latitude = runif(20, min = 34.0, max = 45.0),
    longitude = runif(20, min = -120.0, max = -75.0),
    name = paste("Sample", 1:20),
    color_by = c(rep('type;proportion_infected;comically_small_test_column', 10), rep('type', 10)),
    type = sample(c("Nursery", "Forest", "Urban", "Farm"), 20, replace = TRUE),
    proportion_infected = runif(20),
    comically_small_test_column = sample(c(0.000001898, 0.0000000023678876, 0.0024), 20, replace = TRUE)
  )
  
  # Define which columns can be used for color or size
  split_colorby = strsplit(metadata$color_by, split = ';')
  unlisted_colorby = unlist(split_colorby)
  plot_factors = unique(unlisted_colorby)
  
  # Create human readable titles for columns
  title_key = tools::toTitleCase(gsub('_', ' ', plot_factors))
  names(title_key) <- plot_factors
  
  # Define what columns can be used for size and color
  color_encodable <- plot_factors
  can_be_encoded_by_size <- vapply(metadata[plot_factors], is.numeric, FUN.VALUE = logical(1))
  size_encodable <- plot_factors[can_be_encoded_by_size]
  
  # Define default color and size columns
  color_column <- NA
  size_column <- NA
  if (length(plot_factors) == 1) {
    color_column <- plot_factors
  } else {
    if (length(size_encodable) >= 1) {
      size_column <- size_encodable[1]
      color_column <- plot_factors[plot_factors != size_column][1]
    } else {
      color_column <- plot_factors[1]
    }
  }

  # Create popup HTML
  make_popup <- function(index, max_decimals = 3) {
    values <- sapply(metadata[index, plot_factors], function(x) {
      if (is.numeric(x)) {
        return(signif(x, digits = max_decimals))
      } else {
        return(x)
      }
    })
    paste0(
      "<b>", metadata$name[index], "</b><br>",
      paste(title_key, ":", values, collapse = "<br>")
    )
  }
  metadata$popup <- unlist(lapply(1:nrow(metadata), make_popup))
  
  
  make_palette <- function(column) {
    values <- metadata[[column]]
    if (is.numeric(values)) {
      return(leaflet::colorNumeric(palette = "viridis", domain = values))
    } else {
      return(leaflet::colorFactor(palette = "viridis", domain = values))
    }
  }
  

  # Create map widget with leaflet
  color_palette_func <- make_palette(color_column)
  map_widget <- leaflet::leaflet(data = metadata)
  map_widget <- leaflet::addProviderTiles(map_widget, leaflet::providers$CartoDB.Positron)
  map_widget <- leaflet::addCircleMarkers(
      map_widget,
      lng = metadata$longitude,
      lat = metadata$latitude,
      clusterOptions = leaflet::markerClusterOptions(),
      color = color_palette_func(metadata[[color_column]]),
      fillOpacity = 0.7,
      stroke = TRUE,
      popup = ~popup,
      label = metadata$sample_id)
  
  map_widget <- leaflet::addLegend(
    map_widget,
    "bottomright", 
    pal = color_palette_func, 
    values = metadata[[color_column]],
    title = unname(title_key[color_column]),
    opacity = 1
  )
  

  
  return(map_widget)
}

  