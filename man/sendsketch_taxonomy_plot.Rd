% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotters.R
\name{sendsketch_taxonomy_plot}
\alias{sendsketch_taxonomy_plot}
\title{Make sunburst plot of sendsketch taxonomy}
\usage{
sendsketch_taxonomy_plot(path, interactive = FALSE, ...)
}
\arguments{
\item{path}{The path to one or more folders that contain
pathogensurveillance output or a table in the format of the
\code{\link[=sendsketch_parsed]{sendsketch_parsed()}} output.}

\item{interactive}{Whether or not to produce an interactive
HTML/javascript-based plot or a static one.}

\item{...}{Passed to \code{sendsketch_best_hits}}
}
\description{
Converts classifications of top hits in sendsketch output into an interactive
sunburst plot.
}
\examples{
path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
sendsketch_taxonomy_plot(path)
sendsketch_taxonomy_plot(path, interactive = TRUE)
\dontshow{
## R CMD check: make sure any open connections are closed afterward
closeAllConnections()
}
}
