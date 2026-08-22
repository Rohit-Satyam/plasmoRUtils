test_that("searchApicoTFdb works", {
  df <- searchApicoTFdb(org="pf")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchMidb works", {
  df <- searchMidb(midbSpecies$`Available Species`[196])
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchMiip works", {
  df <- searchMiip(c("PF3D7_0807800","PF3D7_1023900"))
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchIpDb works", {
  df <-  searchIpDb( c("PF3D7_0807800", "PF3D7_1023900"))
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchKipho works", {
  df <-  searchKipho(org = "pf", type = "kinase")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchHP works", {
  vcr::use_cassette("hitpredict",{
    df <- searchHP("Q8I1Q4")
  })
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchRS works", {
    df <-  searchRS(
      geneID = "PfAP2-I",
      org = "Plasmodium falciparum",
      gene_aliases = c(
        "Pf AP2-I",
        "pfap2-i",
        "PF3D7_1007700",
        "Apetala 2 Invasion",
        "AP2-I transcription factor"
      ))
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchPM works", {
    df <- searchPM(geneID = c("PF3D7_0420300","PF3D7_0621000"))
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchTedConsensus works", {
    df <-  searchTedConsensus(c("Q7K6A1","Q8IAP8","C0H4D0","C6KT90","Q8IBJ7"),returnCATHdesc=FALSE)

  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})

test_that("searchPhPl works", {
  vcr::use_cassette("phenoplasm",{
    df <-  searchPhPl(geneID = c("PBANKA_0413500"), org="pb")
  })
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 1)
})