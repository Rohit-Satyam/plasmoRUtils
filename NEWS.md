# plasmoRUtils 1.1.3

In the latest release (1.1.3), an initial static code audit of plasmoRUtils v1.1.2 was conducted with assistance from OpenAI's ChatGPT (GPT-6, October 2026). The tool was used to examine the R source code, identify potential bugs and implementation oversights, and suggest corrective changes.

All AI-generated findings and proposed modifications were subject to review by the package maintainer, who remains responsible for the final implementation, validation, and maintenance of the software.

## Changes in version 1.1.3

### New
- `getTableVeupathdb()` enables you to fetch multiple fields from VEuPathDB given a vector of mixed identifiers from different species.
- `toGeneid()` now uses `getTableVeupathdb()` function to perform desired ID mappings. Unidirectional mapping has been replaced with bidrectional mappings. Organism name and database name is no longer required.

### Major fixes
- `easyTPM` lacked reordering of gene info data frame w.r.t count matrix. This would lead to incorrect TPM normalization if user supplied gene lengths data frame without reordering rows w.r.t count matrix. This scenario is now handled.
- `easytopGO`: The background geneset was being passed but was not used in enrichment analysis. This has been fixed now. A bug in the visualization function `easyGOPlot` was also fixed.

### Documentation and vignette updates
- Revised the tutorials against the version 1.1.3 function signatures, including explicit Fisher/KS input requirements, GO background, API-gated examples and non-ambiguous legacy/probe mappings.
- Updated the README and vignette dependencies and corrected unsupported examples.

### Minor fixes
- `easyGOPlot` legend was being plotted incorrectly for some borderline cases. This has now been fixed.
- Improved `easymakeOrgdb` function. Negative (`NOT`) GO annotations, invalid/unknown GO identifiers, and associations whose gene IDs do not match the GFF are excluded. API access is implemented.
- API access also implemented in `easyTxDbmaker`.
- Fixed temporary file creation and deletion in `easyMCA`.

## Changes in version 1.1.2

- VEuPathDB databases related functions were rewritten to enable usage of API keys for access. These include `getPreconfiguredTable`, `getTable`, `getpairedOrthologs`, `listOrthomcl` and `listVeupathdb`.
- Added new function `getTableAttributes` to list attributes for `customFields` argument in `getTable` function.
- Added caching capability in `getPreconfiguredTableOrthomcl` to reduce time for massive query and subsetting.
- Fixed linking functions across documented help pages.

## Changes in version 1.1.1

- Improved error handling in `searchTedConsensus` function for marginal cases where uniprot IDs throws 404 error.
- Error handling improved in `easyAF2Signal` for marginal cases like `A0A140KXW6`.
- Deprecation of `searchGSC` due to compliance issues. This function has now been replaced with `searchRS` to query Research Square programmatically using its API.
- Removed `ggpubr`, `ggsci` and `glue` to slim down package.
- New function: `getPreconfiguredTableOrthomcl()` is a new function to fetch preconfigured tables from OrthoMCL. It can be used to map old OG IDs to new OG IDs, retrive tables of PFam Architecture of Each Protein, Summary of Pfam domains, Summary of EC Numbers etc.
- Fixed latency in `searchHP` by using `httr2`.
- Revised `searchPM` function by replacing depreciated functions with new functions. 
- Revised `searchPhPl` to handle borderline cases.


## Changes in version 1.1.0

-   `searchGSC` function can now fetch articles for custom years of interest. Title translation to desired language provided. Users can choose how many pages to parse to match number of articles returned upon manual query in Google scholar.
-   `searchHP` : Taxon id is not mandatory but desired to avoid fetching interactions from other organisms having same gene name.
-   Help sections of functions improved.
-   **New functions:** `getpairedOrthologs` , `listOrthomcl`, `listipdb` to get paired orthologs from OrthoMCL 7 and InParanoiDB.
-   **Functions removed:** `searchMT` and `easyPie` were removed due to repeated failure given the database latency.
-   `searchApicoTFdb` now fetches more specialised tables from the database. `searchIpDb` is not *Plasmodium* specific and can accept uniprot IDs of other organisms provided by InParanoiDB.
-   `searchPhPl` was unable to fetch phenotypes for some genes. This has been fixed.

## Changes in version 1.0.0

-   First stable release.

## Changes in version 0.99.1

-   Vignettes created.

## Changes in version 0.99.0

-   First version created.
