#' Pick assemblies to provide context for one or more taxa
#'
#' Pick assemblies to provide context for one or more taxa. For each taxon
#' supplied, representative
#'
#' @param path The path to one or more folders that contain pathogensurveillance
#'   output. The needed files will be automatically found. Alternatively, the
#'   `metadata` and `taxon` parameters can be used.
#' @param metadata `data.frame`/`tibble`: Metadata on assemblies that can be
#'   chosen. This is an alternative to the path `path` input.
#' @param taxon `character` vector: One or more taxa to find assemblies for.
#'   This is an alternative to the `path` input.
#' @param rank `character` vector: The rank corresponding to each value in the
#'   `taxon` argument. Generally not needed but can be used to avoid unexpected
#'   behavior from the same taxon name being used in multiple ranks. Must be
#'   composed of values present in `rank_col`. If not supplied, taxon names
#'   appearing in multiple ranks will be treated as different taxa.
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
#'
#' @return A copy of `metadata` with additional columns `selected`,
#'   `selection_rank`, `selection_taxon`, and  `selection_subtaxon`.
#'
#' @keywords internal
pick_assemblies <- function(
    path = NULL,
    metadata = NULL,
    taxon = NULL,
    rank = NULL,
    subtaxon_count = 5,
    assembly_count = 1,
    rank_col = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'name'),
    criteria = c('is_typical', 'is_type', 'is_refseq', 'is_binomial', 'is_annotated',
                 'assembly_level', 'completeness', '-contamination', 'contig_l50', 'coverage'),
    allow_ambiguous = TRUE, 
    allow_atypical = TRUE
) {
  if (is.null(path)) {
    
  } else {
    metadata <- considered_ref_meta_parsed(path)
  }
  
  # Filter out references with non-standard names
  is_ambiguous <- function(x) {
    ambiguous_words <- c(
      'uncultured',
      'unknown',
      'incertae sedis',
      'sp\\.',
      'cf\\.',
      'endosymbiont',
      'symbiont'
    )
    ambiguous_pattern <- paste0('\\b', ambiguous_words, '\\b', collapse = '|')
    grepl(x, pattern = ambiguous_pattern, ignore.case = TRUE)
  }
  is_latin_binomial <- function(x) {
    grepl(x, pattern = '^[a-zA-Z]+ [a-zA-Z]+($| ).*$') & ! is_ambiguous(x)
  }
  if (only_binomial) {
    assem_data <- assem_data[is_latin_binomial(assem_data$species), ]
  }
  
  # Parse "count" arguments which can be a number or a percentage
  get_count <- function(choices, count) {
    if (grepl(count, pattern = "%$")) {
      prop <- as.numeric(sub(count, pattern = "%", replacement = "")) / 100
      count <- ceiling(choices * prop)
      return(min(c(choices, count)))
    } else {
      count <- as.numeric(count)
      return(min(c(choices, count)))
    }
  }
  
  # Sort references by desirability
  priority <- order(
    decreasing = TRUE,
    assem_data$is_atypical == FALSE,
    assem_data$is_type, # Is type strain
    assem_data$source_database == 'SOURCE_DATABASE_REFSEQ', # Is a RefSeq reference
    is_latin_binomial(assem_data$species), # Has a species epithet
    assem_data$is_annotated,
    factor(assem_data$assembly_level, levels = c("Contig", "Scaffold", "Chromosome", "Complete Genome"), ordered = TRUE),
    assem_data$checkm_completeness,
    -1 * assem_data$checkm_contamination,
    assem_data$contig_l50,
    assem_data$coverage
  )
  assem_data <- assem_data[priority, ]
  
  # Initialize column to hold which level an assembly is selected for
  assem_data$selection_rank <- NA
  assem_data$selection_taxon <- NA
  assem_data$selection_subtaxon <- NA
  
  # Select representatives for each rank
  select_for_rank <- function(assem_data, query_taxa, rank, subrank, count_per_rank, count_per_subrank = 1)  {
    for (tax in query_taxa) {
      # Get assembly indexes for every subtaxon
      tax_to_consider <- (assem_data[[rank]] == tax | assem_data[[subrank]] == tax) & is.na(assem_data$selection_rank)
      subtaxa_found <- unique(assem_data[[subrank]][tax_to_consider])
      subtaxa_found <- subtaxa_found[! subtaxa_found %in% c(assem_data$selection_taxon, assem_data$selection_subtaxon)] # Dont include the data for taxa already chosen
      selected <- lapply(subtaxa_found, function(subtax) {
        which(assem_data[[subrank]] == subtax & is.na(assem_data$selection_rank))
      })
      names(selected) <- subtaxa_found
      
      # Parse count attributes, which can be percentages or integers
      count_per_rank <- get_count(length(selected), count_per_rank)
      count_per_subrank <- get_count(length(selected), count_per_subrank)
      
      # Pick subtaxa with the most assemblies and best mean attributes (based on order in input)
      mean_index <- vapply(selected, mean, FUN.VALUE = numeric(1))
      subtaxa_count <- vapply(selected, length, FUN.VALUE = numeric(1))
      selection_priority <- order(decreasing = TRUE,
                                  is_ambiguous(names(selected)) == FALSE,
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
    return(assem_data)
  }
  assem_data <- select_for_rank(
    assem_data,
    query_taxa = species,
    rank = 'species',
    subrank = 'organism_name',
    count_per_rank = n_ref_strains
  )
  assem_data <- select_for_rank(
    assem_data,
    query_taxa = genera,
    rank = 'genus',
    subrank = 'species',
    count_per_rank = n_ref_species
  )
  assem_data <- select_for_rank(
    assem_data,
    query_taxa = families,
    rank = 'family',
    subrank = 'genus',
    count_per_rank = n_ref_genera
  )
  
  result <- assem_data[! is.na(assem_data$selection_taxon), ]
  
  # Reformat results to the same format as the user-defined metadata
  if (nrow(result) == 0) {
    formatted_result <- data.frame(
      ref_id = character(0),
      ref_name = character(0),
      ref_description = character(0),
      ref_path = character(0),
      ref_ncbi_accession = character(0),
      ref_ncbi_query = character(0),
      ref_ncbi_query_max = character(0),
      ref_primary_usage = character(0),
      ref_contextual_usage = character(0),
      ref_color_by = character(0)
    )
  } else {
    suffix <- paste0(
      ifelse(result$is_type, 'T', ''),
      ifelse(result$source_database == 'SOURCE_DATABASE_REFSEQ', 'R', ''),
      ifelse(result$is_atypical, 'A', '')
    )
    formatted_result <- data.frame(
      ref_id = result$reference_id,
      ref_name = result$organism_name,
      ref_description = paste0(
        result$organism_name, ' ',
        result$accession,
        ifelse(nchar(suffix) > 0, paste0(' ', suffix), '')
      ),
      ref_path = '',
      ref_ncbi_accession = result$accession,
      ref_ncbi_query = '',
      ref_ncbi_query_max = '',
      ref_primary_usage = 'optional',
      ref_contextual_usage = 'optional',
      ref_color_by = ''
    )
  }
  
  # Save to output file
  write.table(result, file = 'raw_results.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(formatted_result, file = out_path, sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(assem_data, file = 'merged_assembly_stats.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
  
  
}




