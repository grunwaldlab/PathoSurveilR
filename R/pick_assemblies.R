#' Pick assemblies to provide context for one or more taxa
#'
#' Pick assemblies to provide context for one or more taxa. For each taxon
#' supplied, representative
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output. The needed files will be automatically found. Alternatively, the
#'   `metadata` and `taxon` parameters can be used.
#' @param subtaxon_count `integer`/`character` vector: The number or percentage
#'   (i.e., a `character` ending with `%`) of taxa to get assemblies for one
#'   rank more specific than the rank of the corresponding value in `taxon`.
#'   This parameter can either be one value for each value in `taxon` or a
#'   single value for all taxa. For example, if `subtaxon_count = c('3', 100%')`
#'   and `taxon = c('Escherichia', 'Lactobacillaceae')` then assemblies for 3
#'   species within Escherichia and all genera within Lactobacillaceae will be
#'   returned.
#' @param assembly_count `integer`/`character` vector: The number or percentage
#'   (i.e., a `character` ending with `%`) of assemblies returned to represent
#'   each subtaxon. This parameter can either be one value for each value in
#'   `taxon` or a single value for all taxa.
#' @param rank_col `character` vector: The names of columns in `metadata`
#'   containing taxon names at each rank, in order from broad to specific (e.g.
#'   `c('family', 'genus', 'species')`). Should include one rank more specific
#'   than any of the taxa defined in `taxon`.
#' @param criteria `character` vector or `list`: The criteria used to decide
#'   which assemblies are the best. If a `character` vector then values are
#'   column names in `metadata` in which higher numbers are better and `TRUE` is
#'   better than `FALSE`. Prefixing the column name with `-` will reverse this
#'   making lower numbers better, etc. If a `list` is supplied, then each item
#'   should be a vector corresponding to rows in `metadata` in which higher
#'   numbers are better and `TRUE` is better than `FALSE`. In both cases,
#'   elements should be sorted in order of decreasing importance.
#' @param allow_ambiguous Allow assemblies with things like 'uncultured',
#'   'unknown', 'sp.', or 'incertea sedis' in taxon names.
#' @param allow_atypical Allow assemblies with numbers or unusual characters in
#'   taxon names.
#' @inheritParams postprocess_table_list
#'
#' @return A copy of `metadata` with additional columns `selected`,
#'   `selection_rank`, `selection_taxon`, and  `selection_subtaxon`.
#'
#' @export
pick_assemblies <- function(
    path = NULL,
    subtaxon_count = 5,
    assembly_count = 1,
    rank_col = c('family', 'genus', 'species', 'organism_name'),
    criteria = c('-is_atypical', 'is_type', 'is_refseq', 'is_binomial', 'is_annotated',
                 'assembly_level', 'checkm_completeness', '-checkm_contamination', '-contig_l50', 'coverage'),
    allow_ambiguous = TRUE, 
    allow_atypical = TRUE,
    simplify = TRUE
) {
  all_assem_data <- considered_ref_meta_parsed(path, simplify = FALSE)
  all_taxa_found <- sendsketch_taxa_present(path, simplify = FALSE)
  
  output <- lapply(names(all_assem_data), function(outdir_path) {
    assem_data <- all_assem_data[[outdir_path]]
    taxa_found <- all_taxa_found[[outdir_path]]
    taxon_id_key <- taxa_found$name
    names(taxon_id_key) <- taxa_found$id
    assem_data$family <- taxon_id_key[assem_data$family_id] 
    
    # Filter out references with non-standard names
    if (! allow_ambiguous) {
      assem_data <- assem_data[is_latin_binomial(assem_data$species), , drop = FALSE]
    }
    
    # Sort references by desirability
    assem_data <- sort_df_by_columns(assem_data, criteria)
    
    # Initialize column to hold which level an assembly is selected for
    assem_data$selection_rank <- NA
    assem_data$selection_taxon <- NA
    assem_data$selection_subtaxon <- NA
    
    # Select representatives for each rank
    taxa_found <- taxa_found[order(match(taxa_found$rank, rank_col), decreasing = TRUE), , drop = FALSE]
    for (i in seq_len(nrow(taxa_found))) {
      tax = taxa_found$name[i]
      rank = taxa_found$rank[i]
      subrank = rank_col[which(rank_col == rank) + 1]
      
      # Get assembly indexes for every subtaxon
      tax_to_consider <- (assem_data[[rank]] == tax | assem_data[[subrank]] == tax) & is.na(assem_data$selection_rank)
      subtaxa_found <- unique(assem_data[[subrank]][tax_to_consider])
      subtaxa_found <- subtaxa_found[! subtaxa_found %in% c(assem_data$selection_taxon, assem_data$selection_subtaxon)] # Dont include the data for taxa already chosen
      selected <- lapply(subtaxa_found, function(subtax) {
        which(assem_data[[subrank]] == subtax & is.na(assem_data$selection_rank))
      })
      names(selected) <- subtaxa_found
      
      # Parse count attributes, which can be percentages or integers
      count_per_rank <- get_count(length(selected), subtaxon_count)
      count_per_subrank <- get_count(length(selected), assembly_count)
      
      # Pick subtaxa with the most assemblies and best mean attributes (based on order in input)
      mean_index <- vapply(selected, mean, FUN.VALUE = numeric(1))
      subtaxa_count <- vapply(selected, length, FUN.VALUE = numeric(1))
      selection_priority <- order(decreasing = TRUE,
                                  is_ambiguous_taxon(names(selected)) == FALSE,
                                  subtaxa_count,
                                  -mean_index
      )
      selected <- selected[selection_priority]
      selected <- selected[seq_len(min(c(count_per_rank, length(selected))))]
      
      # Pick representatives of subtaxa with best attributes (based on order in input)
      selected <- lapply(selected, function(x) {
        x[seq_len(min(c(count_per_subrank, length(x))))]
      })
      
      # Record data on selected assemblies
      selected <- unlist(selected)
      assem_data$selection_rank[selected] <- rank
      assem_data$selection_taxon[selected] <- tax
      assem_data$selection_subtaxon[selected] <- assem_data[[subrank]][selected]
    }
    
    assem_data[! is.na(assem_data$selection_taxon), ]
  })
  names(output) <- names(all_assem_data)
  
  postprocess_table_list(output, simplify = simplify)
}




