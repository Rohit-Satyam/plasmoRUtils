# Changelog

## plasmoRUtils 1.1.1

### Changes in version 1.1.1

- Improved error handling in `searchTedConsensus` function for marginal
  cases where uniprot IDs throws 404 error.
- Error handling improved in `easyAF2Signal` for marginal cases like
  `A0A140KXW6`.
- Deprecation of `searchGSC` due to compliance issues. This function has
  now been replaced with `searchRS` to query Research Square
  programmatically using its API.
- Removed `ggpubr`, `ggsci` and `glue` to slim down package.
- New function:
  [`getPreconfiguredTableOrthomcl()`](https://rohit-satyam.github.io/plasmoRUtils/reference/getPreconfiguredTableOrthomcl.md)
  is a new function to fetch preconfigured tables from OrthoMCL. It can
  be used to map old OG IDs to new OG IDs, retrive tables of PFam
  Architecture of Each Protein, Summary of Pfam domains, Summary of EC
  Numbers etc.
- Fixed latency in `searchHP` by using `httr2`.
- Revised `searchPM` function by replacing depreciated functions with
  new functions.
- Revised `searchPhPl` to handle borderline cases.

### Changes in version 1.1.0

- `searchGSC` function can now fetch articles for custom years of
  interest. Title translation to desired language provided. Users can
  choose how many pages to parse to match number of articles returned
  upon manual query in Google scholar.
- `searchHP` : Taxon id is not mandatory but desired to avoid fetching
  interactions from other organisms having same gene name.
- Help sections of functions improved.
- **New functions:** `getpairedOrthologs` , `listOrthomcl`, `listipdb`
  to get paired orthologs from OrthoMCL 7 and InParanoiDB.
- **Functions removed:** `searchMT` and `easyPie` were removed due to
  repeated failure given the database latency.
- `searchApicoTFdb` now fetches more specialised tables from the
  database. `searchIpDb` is not *Plasmodium* specific and can accept
  uniprot IDs of other organisms provided by InParanoiDB.
- `searchPhPl` was unable to fetch phenotypes for some genes. This has
  been fixed.

### Changes in version 1.0.0

- First stable release.

### Changes in version 0.99.1

- Vignettes created.

### Changes in version 0.99.0

- First version created.
