# These query live services, so failures may be caused by remote access.

# ---- GET functions (VEuPathDB/OrthoMCL API key required) ----
test_that("getTableAttributes works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- getTableAttributes(db = "plasmodb")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("getTable works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- getTable("Plasmodium falciparum 3D7", db = "plasmodb")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("getPreconfiguredTable works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- getPreconfiguredTable("Plasmodium falciparum 3D7",
                              db = "plasmodb", customField = "Y2hInteractions")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("getPreconfiguredTableOrthomcl works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- getPreconfiguredTableOrthomcl("OG3_10277")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)  # The queried group can legitimately have one row.
})

test_that("getpairedOrthologs works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  ids <- listOrthomcl()
  from <- ids$ID[ids$Organism == "Plasmodium falciparum 3D7"]
  to <- ids$ID[ids$Organism == "Toxoplasma gondii ME49"]
  expect_gt(length(from), 0)
  expect_gt(length(to), 0)
  df <- getpairedOrthologs(from = from[1], to = to[1], db = "orthomcl")
  expect_s3_class(df, "data.frame")
})

test_that("toGeneid works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- toGeneid("PF3D7_0420300", from = "ensembl")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

# ---- LIST functions ----
test_that("listOrthomcl works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- listOrthomcl()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("listVeupathdb works", {
  
  skip_if(!nzchar(Sys.getenv("VEUPATHDB_API_KEY")), "VEuPathDB API key missing")
  df <- listVeupathdb()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("listipdb works", {
  
  df <- listipdb()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("pdb2uniprot works", {
  
  df <- pdb2uniprot("4HHB")  # Known haemoglobin PDB with UniProt mappings.
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

# ---- SEARCH functions ----
test_that("searchApicoTFdb works", {
  
  df <- searchApicoTFdb(org = "pf")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("searchMidb works", {
  
  df <- searchMidb(midbSpecies$`Available Species`[196])
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("searchMiip works", {
  
  df <- searchMiip(c("PF3D7_0807800", "PF3D7_1023900"))
  expect_s3_class(df, "data.frame")
})

# test_that("searchIpDb works", {
#   
#   df <- suppressWarnings(searchIpDb("Q8IAR6", idtype = "uniprot"))
#   expect_s3_class(df, "data.frame")
#   expect_gt(nrow(df), 0)
# })

test_that("searchKipho works", {
  
  df <- searchKipho(org = "pf", type = "kinase")
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("searchHP works", {
  
  # Suppress the optional VEuPathDB conversion; exercise the public fetch.
  df <- searchHP("Q8I1Q4", uniprotToGID = FALSE)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("searchRS works", {
  
  df <- searchRS("PfAP2-I", "Plasmodium falciparum",
                 gene_aliases = "PF3D7_1007700", max_pages = 1, sleep_sec = 0)
  expect_s3_class(df, "data.frame")
})

test_that("searchPM works", {
  
  df <- searchPM("PF3D7_0420300", verbose = FALSE)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("searchTedConsensus works", {
  
  df <- searchTedConsensus("Q7K6A1", returnCATHdesc = FALSE)
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

# test_that("searchPhPl works", {
#   
#   df <- searchPhPl("PBANKA_0413500", org = "pb")
#   expect_s3_class(df, "data.frame")
#   expect_gt(nrow(df), 0)
# })
