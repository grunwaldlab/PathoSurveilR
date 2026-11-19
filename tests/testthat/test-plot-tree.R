make_plot_tree_input <- function(tip_counts) {
  trees <- lapply(seq_along(tip_counts), function(i) {
    tip_ids <- paste0("sample-", i, "-", seq_len(tip_counts[[i]]))
    tree <- ape::rtree(tip_counts[[i]])
    tree$tip.label <- tip_ids
    tree
  })
  names(trees) <- paste0("tree-", seq_along(trees))
  tip_ids <- unlist(lapply(trees, `[[`, "tip.label"), use.names = FALSE)
  sample_meta <- data.frame(
    sample_id = tip_ids,
    name = tip_ids,
    color_by = "",
    stringsAsFactors = FALSE
  )
  ref_meta <- data.frame(
    ref_id = "reference-1",
    ref_name = "Reference 1",
    ref_description = "",
    stringsAsFactors = FALSE
  )

  list(tree = trees, sample_meta = sample_meta, ref_meta = ref_meta)
}

test_that("plot_tree height scales to the largest tree", {
  input <- make_plot_tree_input(c(10L, 50L))

  widget <- PathoSurveilR::plot_tree(
    tree = input$tree,
    sample_meta = input$sample_meta,
    ref_meta = input$ref_meta
  )

  expect_equal(widget$height, 1000)
  expect_length(widget$x$trees, 2)
})

test_that("plot_tree height is bounded for small and large trees", {
  small_input <- make_plot_tree_input(5L)
  small_widget <- PathoSurveilR::plot_tree(
    tree = small_input$tree,
    sample_meta = small_input$sample_meta,
    ref_meta = small_input$ref_meta
  )
  expect_equal(small_widget$height, 700)

  large_input <- make_plot_tree_input(100L)
  large_widget <- PathoSurveilR::plot_tree(
    tree = large_input$tree,
    sample_meta = large_input$sample_meta,
    ref_meta = large_input$ref_meta
  )
  expect_equal(large_widget$height, 1600)
})

test_that("plot_tree accepts an explicit widget height", {
  input <- make_plot_tree_input(10L)

  widget <- PathoSurveilR::plot_tree(
    tree = input$tree,
    sample_meta = input$sample_meta,
    ref_meta = input$ref_meta,
    height = "900px"
  )

  expect_identical(widget$height, "900px")
})
