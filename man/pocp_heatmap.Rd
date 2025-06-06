% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotters.R
\name{pocp_heatmap}
\alias{pocp_heatmap}
\title{Plot POCP matrix}
\usage{
pocp_heatmap(
  path,
  combine = FALSE,
  interactive = FALSE,
  height = NULL,
  width = NULL,
  dpi = 100,
  font_size = 8,
  max_label_length = 30
)
}
\arguments{
\item{path}{The path to one or more folders that contain pathogensurveillance
output or paths to tree files.}

\item{combine}{If \code{TRUE}, combine data from all POCP matrices found into a
single plot.}

\item{interactive}{Whether to use an HTML-based interactive format or not
(default: TRUE)}

\item{height}{The height in pixels. If not \code{interactive}, this is divided by
\code{dpi} to convert it to inches.}

\item{width}{The width in pixels. If not \code{interactive}, this is divided by
\code{dpi} to convert it to inches.}

\item{dpi}{How pixels are converted to inches}

\item{font_size}{Size of text used for labels}

\item{max_label_length}{Labels longer than this length will be shortened.}
}
\value{
A list of plots, unless \code{combine} is used, in which case a single
plot is returned.
}
\description{
Plot POCP matrix with dendrogram from data present in the output of a
pathogensurveillance run.
}
\examples{
path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
pocp_heatmap(path)
pocp_heatmap(path, interactive = TRUE)

}
