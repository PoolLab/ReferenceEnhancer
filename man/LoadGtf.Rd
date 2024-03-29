% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadGtf.R
\name{LoadGtf}
\alias{LoadGtf}
\title{LoadGtf}
\usage{
LoadGtf(unoptimized_annotation_path)
}
\arguments{
\item{unoptimized_annotation_path}{Path to the unoptimized genome annotion GTF file.}
}
\value{
Resulting object contains the genome annotation entries from the genome annotation GTF file.
}
\description{
Use to import the Ensembl/10x Genomics default genome annotation
or other desired genome annotation file in GTF format for optimization for scRNA-seq
analysis. Note: This file can be downloaded from 10x Genomics provided reference
transcriptome "gene" folder at
"https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"
or Ensembl.org if wish to customize more.
}
\examples{
LoadGtf(unoptimized_annotation_path = "test_genes.gtf")
}
