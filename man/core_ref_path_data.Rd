% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_data_finders.R
\name{core_ref_path_data}
\alias{core_ref_path_data}
\title{Find the core gene analysis reference path data}
\usage{
core_ref_path_data(path)
}
\arguments{
\item{path}{The path to one or more folders that contain
pathogensurveillance output.}
}
\value{
\code{tibble} with \code{report_group_id} and \code{path} columns
}
\description{
Return a table with the file path to the TSV with the list of references used
in the core gene analysis for a given pathogensurveillance output folder.
}
\examples{
path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
core_ref_path_data(path)

}
\seealso{
Other path tables: 
\code{\link{busco_ref_path_data}()},
\code{\link{busco_tree_path_data}()},
\code{\link{considered_ref_meta_path_data}()},
\code{\link{core_tree_path_data}()},
\code{\link{estimated_ani_matrix_path_data}()},
\code{\link{pocp_matrix_path_data}()},
\code{\link{ref_meta_path_data}()},
\code{\link{run_info_path_data}()},
\code{\link{sample_meta_path_data}()},
\code{\link{selected_ref_meta_path_data}()},
\code{\link{sendsketch_path_data}()},
\code{\link{software_version_path_data}()},
\code{\link{status_message_path_data}()},
\code{\link{variant_align_path_data}()},
\code{\link{variant_ref_path_data}()},
\code{\link{variant_tree_path_data}()}
}
\concept{path tables}
