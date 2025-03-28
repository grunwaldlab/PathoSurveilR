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
  
  ##### CODE TO ADAPT #####
  
  library(RcppSimdJson)
  
  # Set random number generator seed
  set.seed(1)
  
  # Parse taxonomy inputs
  args <- commandArgs(trailingOnly = TRUE)
  args <- c(
    '/home/fosterz/projects/pathogensurveillance/work/33/7afb691e09b32a756b9073ce718108/LF1_families.txt',
    '/home/fosterz/projects/pathogensurveillance/work/33/7afb691e09b32a756b9073ce718108/LF1_genera.txt',
    '/home/fosterz/projects/pathogensurveillance/work/33/7afb691e09b32a756b9073ce718108/LF1_species.txt',
    '3', '3', '3', 'false', 'LF1.tsv',
    '/home/fosterz/projects/pathogensurveillance/work/33/7afb691e09b32a756b9073ce718108/Lactobacillaceae.json',
    '/home/fosterz/projects/pathogensurveillance/work/33/7afb691e09b32a756b9073ce718108/Streptomycetaceae.json'
  )
  args <- as.list(args)
  families <- readLines(args[[1]])
  genera <- readLines(args[[2]])
  species <- readLines(args[[3]])
  n_ref_strains <- as.numeric(args[[4]])
  n_ref_species <- as.numeric(args[[5]])
  n_ref_genera <- as.numeric(args[[6]])
  only_binomial <- as.logical(args[[7]])
  out_path <- args[[8]]
  
  # Parse input JSONs
  if (length(args) < 9) {
    stop('No family-level reference metadata files supplied. Check input data.')
  }
  json_paths <- unlist(args[9:length(args)])
  json_families <- gsub(basename(json_paths), pattern = '.json', replacement = '', fixed = TRUE)
  json_data <- lapply(seq_along(json_paths), function(index) {
    parsed_json <- RcppSimdJson::fparse(readLines(json_paths[index]), always_list = TRUE)
    output <- do.call(rbind, lapply(parsed_json, function(assem_data) {
      attributes <- assem_data$assembly_info$biosample$attributes
      hosts <- paste0(attributes$value[attributes$name == 'host'], collapse = ';')
      data_parts <- list(
        accession = assem_data$accession,
        assembly_level = assem_data$assembly_info$assembly_level,
        assembly_status = assem_data$assembly_info$assembly_status,
        assembly_type = assem_data$assembly_info$assembly_type,
        hosts = ifelse(hosts == '', NA_character_, hosts),
        organism_name = assem_data$organism$organism_name,
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
      output$family <- rep(json_families[index], nrow(output))
    }
    return(output)
  })
  assem_data <- do.call(rbind, json_data)
  if (is.null(assem_data)) {
    assem_data <- data.frame(
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
    )
  }
  
  # Add column for modified ID
  modified_id <- gsub(assem_data$accession, pattern = '[\\/:*?"<>| .]', replacement = '_')
  assem_data <- cbind(reference_id = modified_id, assem_data)
  rownames(assem_data) <- NULL
  
  # Add taxon info columns
  assem_data$organism_name <- gsub(assem_data$organism_name, pattern = '[', replacement = '', fixed = TRUE)
  assem_data$organism_name <- gsub(assem_data$organism_name, pattern = ']', replacement = '', fixed = TRUE)
  assem_data$species <- gsub(assem_data$organism_name, pattern = '([a-zA-Z0-9.]+) ([a-zA-Z0-9.]+) (.*)', replacement = '\\1 \\2')
  assem_data$genus <- gsub(assem_data$organism_name, pattern = '([a-zA-Z0-9.]+) (.*)', replacement = '\\1')
  
  # Filter out references with non-standard names
  is_ambiguous <- function(x) {
    ambiguous_words <- c(
      'uncultured',
      'unknown',
      'incertae sedis',
      'sp.'
    )
    vapply(x, FUN.VALUE = logical(1), function(one) {
      any(vapply(ambiguous_words, FUN.VALUE = logical(1), function(word) {
        grepl(one, pattern = word, ignore.case = TRUE)
      }))
    })
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
      
      # Pick subtaxa with the most assemblies and best mean attributes (based on order in input)
      mean_index <- vapply(selected, mean, FUN.VALUE = numeric(1))
      subtaxa_count <- vapply(selected, length, FUN.VALUE = numeric(1))
      selection_priority <- order(
        is_ambiguous(names(selected)),
        -subtaxa_count,
        mean_index
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
    formatted_result <- data.frame(
      ref_id = result$reference_id,
      ref_name = result$organism_name,
      ref_description = paste0(
        result$species, ' (',
        result$accession,
        ifelse(result$is_type, '; T', ''),
        ifelse(result$source_database == 'SOURCE_DATABASE_REFSEQ', '; RS', ''),
        ifelse(result$is_atypical, '; A', ''),
        # ifelse(is.na(result$hosts), '', paste0('; Host: ', result$hosts)),
        ')'
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




