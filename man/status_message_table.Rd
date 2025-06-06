% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table_printers.R
\name{status_message_table}
\alias{status_message_table}
\title{Get table of pipeline status data}
\usage{
status_message_table(input, summarize_by = NULL, interactive = FALSE, ...)
}
\arguments{
\item{input}{The path to one or more folders that contain
pathogensurveillance output.}

\item{summarize_by}{What variable to summarize results by, if any. Can be one
of \code{'sample'}, \code{'message'}, \code{'workflow'}, or \code{NULL}. By default, all values
in the message data is shown on its own row.}

\item{interactive}{Whether to produce interactive tables
(TRUE) or static tables (FALSE). Defaults to TRUE if the environment
supports HTML output, otherwise FALSE. Interactive tables offer enhanced
browsing capabilities, while static tables are best for printed pdf
reports.}

\item{...}{Passed to \code{DT::datatable} for interactive output.}
}
\value{
A table with details for errors, warnings, notes.
}
\description{
Return a formatted interactive table with the data on the issues encountered
by the pipeline, one row for each issue. The contents of all status message
files found in the given paths will be combined.
}
\examples{
path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
status_message_table(path)
status_message_table(path, interactive = TRUE)
status_message_table(path, summarize_by = 'sample')
status_message_table(path, summarize_by = 'workflow')
status_message_table(path, summarize_by = 'message')

}
