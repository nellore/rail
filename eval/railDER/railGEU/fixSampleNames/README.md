This directory has the scripts for matching the sample names from the GEUVADIS
`ballgown` object and the Rail-RNA manifest files. The same script creates the
`GRanges` object with the ERs (`regions.Rdata`) and the coverage matrix
(`coverageMatrix.Rdata`). The latter is a large file, so it's not version
controlled. You can find it via Figshare at [dx.doi.org/10.6084/m9.figshare.1405638](http://dx.doi.org/10.6084/m9.figshare.1405638).

`create_regions_suppTable.R` takes `regions.Rdata` and creates a CSV file with the same information.
