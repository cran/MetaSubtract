\name{alternative_headers}
\alias{alternative_headers}
\title{
A file with alternative header names in the input files
}
\description{
A file with alternative header names in the input files, which is used to 
internally convert the header names from your file(s) to the ones used by
the script. The file has two columns. The first one contains the variable 
names as used in the script, the second one contains the header names as 
specified in the input files. A default file can be found in the subdirectory 
'extdata' of the installation directory of MetaSubtract. You can add header 
names of your files to this file to make the internal conversion.

The variable names in the script are the following:

MARKER       marker name

EFFECTALLELE effect allele

OTHERALLELE  other allele

BETA         effect size

EAF          effect allele frequency

N            sample size

P            p-value

SE           standard error of BETA

Z            z-score

NSTUDIES     number of studies

DIRECTION    direction of effects

LP           -log10(p-value)

QHET         Q-value for heterogeneity

QHETP        p-value of QHET

I2HET        I2 value for heterogeneity

}
