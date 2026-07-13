# Create Org.db package quickly

A convenience function to make OrgDb packages quickly to be used with
other GO enrichment packages such as ClusterProfiler.

## Usage

``` r
easymakeOrgdb(
  gff =
    "https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff",
  gaf =
    "https://plasmodb.org/common/downloads/release-68/Pfalciparum3D7/gaf/PlasmoDB-68_Pfalciparum3D7_Curated_GO.gaf.gz",
  out.dir = ".",
  taxid = 36329,
  genus = "Plasmodium",
  sp = "falciparum3D7",
  version = 0.1,
  verbose = FALSE,
  maintainer = "John doe <johndoe@gmail.com>"
)
```

## Arguments

- gff:

  A link to GFF file from VEuPathDB or path to GFF file if the file is
  present locally.

- gaf:

  Gene Ontology file obtained from VEuPathDB.

- out.dir:

  Output directory where the package will be saved.

- taxid:

  Taxonomy ID of the organism. This can be obtained from
  https://www.ncbi.nlm.nih.gov/taxonomy.

- genus:

  Genus of the organism. This will be used to construct the name of the
  package.

- sp:

  Species with or without strain information.

- version:

  Version of the package if you choose to maintain and share the
  package.

- verbose:

  Display the messages when running makeOrgPackage function.

- maintainer:

  Email Id of the package builder. Default "John doe
  <johndoe@gmail.com>"

## Value

A tar.gz file that can be installed as a package and can be used with GO
enrichment tools such as ClusterProfiler.

## Examples

``` r
if (FALSE) { # \dontrun{
 easymakeOrgdb(
  gff="https://plasmodb.org/common/downloads/release-68/PbergheiANKA/gff/data/PlasmoDB-68_PbergheiANKA.gff",
  gaf="https://plasmodb.org/common/downloads/release-68/PbergheiANKA/gaf/PlasmoDB-68_PbergheiANKA_Curated_GO.gaf.gz",
  out.dir=".", taxid=5823,genus="Plasmodium",
  sp="bergheiANKA",
  version=0.1,
  verbose = FALSE,
  maintainer="John doe <johndoe@gmail.com>")


} # }
```
