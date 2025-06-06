% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotters.R
\name{plot_tree}
\alias{plot_tree}
\title{Plot generic phylogeny}
\usage{
plot_tree(
  path = NULL,
  target = NULL,
  tree = NULL,
  sample_meta = NULL,
  ref_meta = NULL,
  collapse_by_tax = NULL,
  interactive = FALSE
)
}
\arguments{
\item{path}{The path to one or more folders that contain pathogensurveillance
output or paths to tree files.}

\item{target}{The names of one or more output types to search for.}

\item{tree}{One or more trees to plot. Only used if \code{path} is not used.}

\item{sample_meta}{The path to sample metadata or the data itself in the form
of a table. Only used if \code{path} is not used.}

\item{ref_meta}{The path to reference metadata or the data itself in the form
of a table. Only used if \code{path} is not used.}

\item{collapse_by_tax}{A table with taxonomic data encoded as one column per
rank and a \code{sample_id} column with the same values present in
\code{sample_meta}. Used to combine the trees together based on the smallest shared
taxonomic tree. Only used if \code{path} is not used.}

\item{interactive}{Whether to use an HTML-based interactive format or not
(default: TRUE)}
}
\value{
A list of plots, unless \code{collapse_by_tax} is used, in which case a
single plot is returned.
}
\description{
Plot phylogenies present in the output of a pathogensurveillance run.
}
