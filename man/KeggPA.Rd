\name{keggPA}
\alias{keggPA}
\title{Kegg Human Pathway Analysis
Based on Fisher's on tail exact test, overrepresentated pathways are found.}
\usage{
keggPA(hgnc, project_name, kegg.list = NULL)
}
\arguments{
  \item{hgnc}{a vector of Human Gene symbol IDs or entrez
  gene IDs in format of "hsa:XXXX"}

  \item{project_name}{a character stating name of the
  project used in writing the results file in the
  directory}

  \item{kegg.list}{is a dataframe in format of Kegg Pathway
  ID, entrez IDs, Gene Symbol IDs and Name of Pathway
  obtained from kegg directory. By default December 2013
  version will be used}
}
\value{
A dataframe of pathway IDs, pathway names, mapped genes,
p-values and FDRs.
}
\description{
Kegg Human Pathway Analysis Based on Fisher's on tail exact
test, overrepresentated pathways are found.
}
\details{
Based on assumptions of hypergeometric test, we omit all
the genes that are not in kegg database. The hypergeometric
p-value is calcualted based on total genes in annotated in
kegg database.
}
\examples{
keggPA(c("BRCA1","BRCA2","TP53","XRCC5"),project_name="Example_please_delete")
}
\author{
Mehran Karimzadeh
}

