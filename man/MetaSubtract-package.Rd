\name{MetaSubtract-package}
\alias{MetaSubtract-package}
\alias{MetaSubtract}
\docType{package}
\title{
Subtracting Summary Statistics of One or more Cohorts from Meta-GWAS Results
}
\description{
  If results from a meta-GWAS are used for validation in one of the cohorts that was included in the meta-analysis, this will yield biased (i.e. too optimistic) results. The validation cohort needs to be independent from the meta-GWAS results. MetaSubtract will subtract the results of the respective cohort from the meta-GWAS results analytically without having to redo the meta-GWAS analysis using the leave-one-out methodology. It can handle different meta-analyses methods and takes into account if single or double genomic control correction was applied to the original meta-analysis. It can be used for whole GWAS, but also for a limited set of SNPs or other genetic markers.
}
\details{
The core of this package is the function \code{\link{meta.subtract}}.
}
\author{
Ilja M. Nolte

Maintainer: Ilja M. Nolte <i.m.nolte@umcg.nl>
}
\references{
Nolte et al. (2017). Missing heritability: is the gap closing? An analysis of 32 complex traits in the Lifelines Cohort Study. Eur J Hum Genet. 2017;25:877-885.
}
\seealso{ 
 \code{\link{meta.subtract}}.
}
\keyword{ package }
