test_that("variant MSN plotting rejects a static request", {
  expect_error(
    PathoSurveilR:::variant_msn_plot("unused-path", interactive = FALSE),
    "Static MSN plots are not supported"
  )
})

make_variant_msn_fixture <- function(path, fasta_files = list()) {
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
      "    },",
      '    "sample_metadata": {',
      '      "name": "Sample metadata",',
      '      "parsers": {"r": "PathoSurveilR:::parse_tsv"},',
      '      "combiners": {"r": "PathoSurveilR:::combine_tables"},',
      '      "sources": [{"path": "metadata", "pattern": "^sample_metadata\\\\.tsv$"}]',
      "    },",
      '    "reference_metadata": {',
      '      "name": "Reference metadata",',
      '      "parsers": {"r": "PathoSurveilR:::parse_tsv"},',
      '      "combiners": {"r": "PathoSurveilR:::combine_tables"},',
      '      "sources": [{"path": "metadata", "pattern": "^reference_metadata\\\\.tsv$"}]',
      "    }",
      "  }",
      "}"
    ),
    file.path(path, ".pathogensurveillance_output.json")
  )
  writeLines(
    c(
      "sample_id\tname\tgroup",
      "sample-a\tSample A\tcontrol",
      "sample-b\tSample B\tcase",
      "sample-c\tSample C\tcase"
    ),
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

test_that("variant MSN widgets include named FASTA datasets and matching metadata", {
  fasta <- c(
    "all--ref-a.fasta" = ">ref-a\nACGT\n>sample-a\nACGA\n>sample-b\nACGG",
    "all--ref-b.fasta" = ">ref-b\nTTTT\n>sample-a\nTTTA\n>sample-c\nTTTG"
  )
  path <- make_variant_msn_fixture(tempfile("variant-msn-"), as.list(fasta))

  widget <- variant_msn_plot(path, interactive = TRUE)

  expect_s3_class(widget, "msnwidget")
  expect_null(widget$width)
  expect_null(widget$height)
  expect_identical(widget$sizingPolicy$defaultWidth, "100%")
  expect_identical(widget$sizingPolicy$defaultHeight, "70vh")
  expect_false(widget$sizingPolicy$knitr$figure)
  expect_length(widget$x$datasets, 2)
  expect_identical(
    vapply(widget$x$datasets, `[[`, character(1), "name"),
    c(
      "Reference A (ref-a)",
      "Reference B (ref-b)"
    )
  )
  for (dataset in widget$x$datasets) {
    fasta_ids <- sub(
      "^>",
      "",
      grep("^>", strsplit(dataset$fasta, "\n", fixed = TRUE)[[1]], value = TRUE)
    )
    metadata_ids <- vapply(dataset$metadata[[1]]$data, `[[`, character(1), "id")
    expect_false(any(fasta_ids %in% c("ref-a", "ref-b")))
    expect_setequal(fasta_ids, metadata_ids)
  }
})

test_that("variant MSN widget height can be overridden", {
  path <- make_variant_msn_fixture(
    tempfile("variant-msn-height-"),
    list("all--ref-a.fasta" = ">ref-a\nACGT\n>sample-a\nACGA\n>sample-b\nACGG")
  )

  widget <- variant_msn_plot(path, height = "500px")

  expect_identical(widget$sizingPolicy$defaultHeight, "500px")
})

test_that("variant MSN plotting handles a single reference and insufficient input", {
  single <- make_variant_msn_fixture(
    tempfile("variant-msn-single-"),
    list("all--ref-a.fasta" = ">ref-a\nACGT\n>sample-a\nACGA\n>sample-b\nACGG")
  )
  single_widget <- variant_msn_plot(single)
  expect_s3_class(single_widget, "msnwidget")
  expect_length(single_widget$x$datasets, 1)

  too_few <- make_variant_msn_fixture(
    tempfile("variant-msn-few-"),
    list("all--ref-a.fasta" = ">ref-a\nACGT\n>sample-a\nACGA")
  )
  expect_null(variant_msn_plot(too_few))

  no_alignments <- make_variant_msn_fixture(tempfile("variant-msn-empty-"))
  expect_null(variant_msn_plot(no_alignments))
})

test_that("samples without metadata remain in the FASTA dataset", {
  path <- make_variant_msn_fixture(
    tempfile("variant-msn-missing-metadata-"),
    list(
      "all--ref-a.fasta" =
        ">ref-a\nACGT\n>sample-a\nACGA\n>sample-missing\nACGG"
    )
  )

  widget <- variant_msn_plot(path)
  dataset <- widget$x$datasets[[1]]
  fasta_ids <- sub(
    "^>",
    "",
    grep("^>", strsplit(dataset$fasta, "\n", fixed = TRUE)[[1]], value = TRUE)
  )
  metadata <- dataset$metadata[[1]]$data

  expect_identical(fasta_ids, c("sample-a", "sample-missing"))
  expect_identical(vapply(metadata, `[[`, character(1), "id"), fasta_ids)
  expect_true(is.na(metadata[[2]]$name))
})
