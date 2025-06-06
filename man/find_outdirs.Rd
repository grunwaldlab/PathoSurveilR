% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_finders.R
\name{find_outdirs}
\alias{find_outdirs}
\title{Find directories containing pathogensurveillance output}
\usage{
find_outdirs(path = ".", allow_nested = FALSE)
}
\arguments{
\item{path}{\code{character}. The starting paths to search. Defaults to current
directory.}

\item{allow_nested}{\code{logical}. If \code{FALSE}, don't search subdirectories of
directories that already contain the target file. If \code{TRUE}, continue
searching all subdirectories.}
}
\value{
\code{character} vector of paths to directories containing the target
file.
}
\description{
Recursively searches a directory to find the first directory (starting with
the input directory) that contains a file called
'.pathogensurveillance_output.json'.
}
\examples{
\dontrun{
# Search current directory and subdirectories (non-nested)
find_output_dirs()

# Search a specific path with nested searching allowed
find_output_dirs("~/projects", allow_nested = TRUE)
}

}
