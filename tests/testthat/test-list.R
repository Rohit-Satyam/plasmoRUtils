test_that("listveupathdb works", {
    df <- listipdb()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("listmca works", {

    df <- listMCA()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("listOrthomcl works", {
    df <- listOrthomcl()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("listVeupathdb works", {
    df <- listVeupathdb()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("pdb2uniprot works", {
    df <- pdb2uniprot("9FIA")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})