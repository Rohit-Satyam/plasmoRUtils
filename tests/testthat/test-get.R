test_that("getpairedOrthologs works", {
  ids <- listOrthomcl()
  df <- getpairedOrthologs(
    from=ids %>% dplyr::filter(Organism=="Plasmodium falciparum 3D7") %>% dplyr::pull(ID), 
    to=ids %>% dplyr::filter(Organism=="Toxoplasma gondii ME49") %>% dplyr::pull(ID),
    db="orthomcl", 
    transform = TRUE)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("getTable works", {
  
    df <- getTable(org="Plasmodium falciparum 3D7", db="plasmodb")

  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("getPreconfiguredTable works", {

    df <-getPreconfiguredTable(org = "Plasmodium falciparum 3D7",db = "plasmodb",customField = "Y2hInteractions")

  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("getPreconfiguredTableOrthomcl works", {

    df <-getPreconfiguredTableOrthomcl("OG3_10277")

  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("toGeneid works", {

    df <- toGeneid(
      c("PF3D7_0420300", "PF3D7_0621000"),
      from="ensembl")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})