test_that("long-format ANI CSV output parses into a symmetric matrix", {
  pairwise <- data.frame(
    query_name = c("genome-a", "genome-a", "genome-b", "genome-c"),
    match_name = c("genome-b", "genome-c", "genome-c", "genome-d"),
    average_containment_ani = c(0.95, 0.90, 0.80, NA_real_)
  )
  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  utils::write.csv(pairwise, csv_file, row.names = FALSE)

  parser <- eval(parse(text = "PathoSurveilR:::parse_ani_csv"))
  parsed <- parser(csv_file)

  expect_identical(dim(parsed), c(4L, 4L))
  expect_identical(
    rownames(parsed),
    c("genome-a", "genome-b", "genome-c", "genome-d")
  )
  expect_equal(parsed["genome-a", "genome-b"], 0.95)
  expect_equal(parsed["genome-b", "genome-a"], 0.95)
  expect_equal(unname(diag(parsed)), rep(1, 4))
  expect_true(is.na(parsed["genome-c", "genome-d"]))
  expect_true(is.na(parsed["genome-d", "genome-c"]))

  widget <- PathoSurveilR:::make_heatmap(
    parsed,
    data.frame(ref_id = "genome-b", ref_name = "Reference B"),
    data.frame(sample_id = "genome-a", name = "Sample A"),
    interactive = TRUE,
    metric = "ani"
  )
  expect_equal(widget$x$matrix[[1]][[2]], 95)
  expect_identical(widget$x$labels, rownames(parsed))
  expect_null(widget$x$order)
  expect_identical(widget$x$full_labels[1], "Sample A")
  expect_identical(widget$x$metric, "ani")
  expect_identical(widget$x$metric_label, "Estimated ANI similarity (%)")
  expect_true(widget$sizingPolicy$viewer$fill)
  expect_true(widget$sizingPolicy$browser$fill)
})

test_that("ANI CSV parsing reports absent required fields", {
  csv_file <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_file), add = TRUE)
  utils::write.csv(data.frame(query_name = "a"), csv_file, row.names = FALSE)

  expect_error(
    PathoSurveilR:::parse_ani_csv(csv_file),
    "missing required columns: match_name, average_containment_ani"
  )
})
