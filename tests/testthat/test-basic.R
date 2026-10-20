# A fast test that does not require internet access or an API key.
test_that("easyTPM returns correctly scaled expression", {
  counts <- matrix(c(10, 10), ncol = 1,
                   dimnames = list(c("gene1", "gene2"), "sample1"))
  lengths <- data.frame(GeneID = c("gene2", "gene1"),
                        Length = c(2000, 1000))

  result <- easyTPM(counts, lengths)

  expect_s3_class(result, "data.frame")
  expect_equal(sum(result$sample1), 1e6, tolerance = 1e-7)
  expect_gt(result["gene1", "sample1"], result["gene2", "sample1"])
})
