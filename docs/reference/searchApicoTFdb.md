# Fetch data from ApicoTFdb

This function provides ability to query your gene IDs to
[ApicoTFDb](https://bioinfo.icgeb.res.in/PtDB/index.html).

## Usage

``` r
searchApicoTFdb(org = "pf", fetch = "all")
```

## Arguments

- org:

  Abbreviation of organism of interest.

- fetch:

  Describe the tables to be fetched. Default: "all" will fetch all TFs.
  To fetch TRs,CRRs,RNA-regs or Experimentally verified TF, use
  "trs","crrs","rnaregs" and "exptfs" respectively. When using other
  than "all", org argument will be ignored.

  - pb: *Plasmodium berghii*

  - pv: *Plasmodium vivax*

  - pf: *Plasmodium falciparum*

  - pk: *Plasmodium knowlesi*

  - py: *Plasmodium yoelii*

  - pc: *Plasmodium chabaudi*

  &nbsp;

  - tg49: *Toxoplasma Gondii* ME49

  - tg89: *Toxoplasma Gondii* P89

  - cp: *Cryptosporidium parvum*

  - em: *Eimeria maxima*

  - bb: *Babesia bovis*

  - et: *Eimeria tenella*

  - nu: *Neurospora caninum*

  - cy: *Cyclospora cayetanensis*

## Value

df This function returns a dataframe of transcription regulators for
organism of interest from ApicoTFDb.

## Examples

``` r
if (FALSE) { # \dontrun{
test <- searchApicoTFdb(org = "pf")
} # }
```
