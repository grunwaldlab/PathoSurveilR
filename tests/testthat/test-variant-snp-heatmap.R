make_variant_snp_fixture <- function(path, fasta_files = list()) {
  dir.create(path, recursive = TRUE)
  dir.create(file.path(path, "variants"))
  dir.create(file.path(path, "metadata"))
  writeLines(
    c(
      "{",
      '  "outputs": {',
      '    "variants_fasta": {',
      '      "name": "Variant sequences",',
      '      "parsers": {"r": "PathoSurveilR:::parse_fasta"},',
      '      "sources": [{"path": "variants", "pattern": "^(.+?)--(.+?)\\\\.fasta$", "capture_groups": ["report_id", "reference_id"]}]',
      '    },',
      '    "sample_metadata": {',
      '      "name": "Sample metadata",',
      '      "parsers": {"r": "PathoSurveilR:::parse_tsv"},',
      '      "combiners": {"r": "PathoSurveilR:::combine_tables"},',
      '      "sources": [{"path": "metadata", "pattern": "^sample_metadata\\\\.tsv$"}]',
      '    },',
      '    "reference_metadata": {',
      '      "name": "Reference metadata",',
      '      "parsers": {"r": "PathoSurveilR:::parse_tsv"},',
      '      "combiners": {"r": "PathoSurveilR:::combine_tables"},',
      '      "sources": [{"path": "metadata", "pattern": "^reference_metadata\\\\.tsv$"}]',
      '    }',
      '  }',
      '}'
    ),
    file.path(path, ".pathogensurveillance_output.json")
  )
  writeLines(
    c("sample_id\tname", "sample-a\tSample A", "sample-b\tSample B", "sample-c\tSample C"),
    file.path(path, "metadata", "sample_metadata.tsv")
  )
  writeLines(
    c("ref_id\tref_name", "ref-a\tReference A", "ref-b\tReference B"),
    file.path(path, "metadata", "reference_metadata.tsv")
  )
  for (name in names(fasta_files)) {
    writeLines(fasta_files[[name]], file.path(path, "variants", name))
  }
  path
}

test_that("SNP distances count jointly callable nucleotide differences", {
  alignment <- matrix(
    c(
      "A", "C", "G", "T",
      "-", "A", "R", "A",
      "N", "T", "-", "C",
      "A", "C", "G", "T"
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("sample-a", "sample-b", "sample-c", "sample-d"), NULL)
  )

  distance <- PathoSurveilR:::snp_distance_matrix(alignment)

  expect_equal(distance["sample-a", "sample-b"], 2)
  expect_equal(distance["sample-b", "sample-a"], 2)
  expect_equal(distance["sample-a", "sample-c"], 2)
  expect_equal(distance["sample-b", "sample-c"], 2)
  expect_equal(distance["sample-a", "sample-d"], 0)
  expect_identical(unname(diag(distance)), c(0, 0, 0, 0))
})

test_that("SNP distance heatmap combines alignments and preserves unobserved pairs", {
  path <- make_variant_snp_fixture(
    tempfile("variant-snp-"),
    list(
      "report--ref-a.fasta" = ">ref-a\nAAAA\n>sample-a\nAACA\n>sample-b\nAATA",
      "report--ref-b.fasta" = ">ref-b\nTTTT\n>sample-c\nTTCT\n>sample-d\nTTAT"
    )
  )

  widget <- variant_snp_heatmap(path, interactive = TRUE, height = "70vh")
  matrix <- do.call(rbind, lapply(widget$x$matrix, function(row) {
    vapply(row, function(value) {
      if (is.null(value)) NA_real_ else as.numeric(value)
    }, numeric(1))
  }))
  dimnames(matrix) <- list(widget$x$labels, widget$x$labels)

  expect_identical(widget$x$metric, "snp_distance")
  expect_identical(widget$x$metric_label, "SNP differences (count)")
  expect_identical(widget$sizingPolicy$defaultHeight, "70vh")
  expect_equal(matrix["sample-a", "sample-b"], 1)
  expect_equal(matrix["sample-c", "sample-d"], 1)
  expect_true(is.na(matrix["sample-a", "sample-c"]))
  expect_identical(widget$x$full_labels[match("ref-a", widget$x$labels)], "Reference A")
})

test_that("SNP distances mark pairs with no jointly callable bases missing", {
  alignment <- matrix(
    c("-", "N", "A", "R"),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("sample-a", "sample-b"), NULL)
  )

  expect_true(is.na(PathoSurveilR:::snp_distance_matrix(alignment)[1, 2]))
})

test_that("SNP distance heatmap supports static rendering", {
  matrix <- matrix(
    c(0, 2, 2, 0),
    nrow = 2,
    dimnames = list(c("sample-a", "sample-b"), c("sample-a", "sample-b"))
  )
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  }, add = TRUE)

  plot <- PathoSurveilR:::make_heatmap(
    input_matrix = matrix,
    sample_data = NULL,
    ref_data = NULL,
    metric = "snp_distance"
  )

  expect_s3_class(plot, "pheatmap")
})

test_that("SNP heatmap warns if IDs overlap across reference alignments", {
  path <- make_variant_snp_fixture(
    tempfile("variant-snp-overlap-"),
    list(
      "report--ref-a.fasta" = ">ref-a\nAAAA\n>sample-a\nAACA\n>sample-b\nAATA",
      "report--ref-b.fasta" = ">ref-b\nTTTT\n>sample-a\nTTCT\n>sample-c\nTTAT"
    )
  )

  expect_warning(
    variant_snp_heatmap(path, interactive = TRUE),
    "Sequence IDs occur in multiple SNP alignments"
  )
})

test_that("SNP heatmap returns NULL without usable alignments", {
  empty <- make_variant_snp_fixture(tempfile("variant-snp-empty-"))
  one_sequence <- make_variant_snp_fixture(
    tempfile("variant-snp-one-"),
    list("report--ref-a.fasta" = ">ref-a\nAAAA")
  )

  expect_null(variant_snp_heatmap(empty, interactive = TRUE))
  expect_null(variant_snp_heatmap(one_sequence, interactive = TRUE))
})
