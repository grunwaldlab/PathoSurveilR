% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_finders.R
\name{software_version_path}
\alias{software_version_path}
\title{Find the software version file path}
\usage{
software_version_path(path)
}
\arguments{
\item{path}{The path to one or more folders that contain pathogensurveillance output.}
}
\value{
character vector of length 1
}
\description{
Return the file path to the YAML file with the versions of software used for
a given pathogensurveillance output folder.
}
\examples{
path <- system.file('extdata/ps_output', package = 'PathoSurveilR')
software_version_path(path)

}
\seealso{
Other path finders: 
\code{\link{busco_ref_path}()},
\code{\link{busco_tree_path}()},
\code{\link{considered_ref_meta_path}()},
\code{\link{core_ref_path}()},
\code{\link{core_tree_path}()},
\code{\link{estimated_ani_matrix_path}()},
\code{\link{pocp_matrix_path}()},
\code{\link{ref_meta_path}()},
\code{\link{run_info_path}()},
\code{\link{sample_meta_path}()},
\code{\link{selected_ref_meta_path}()},
\code{\link{sendsketch_path}()},
\code{\link{status_message_path}()},
\code{\link{variant_align_path}()},
\code{\link{variant_ref_path}()},
\code{\link{variant_tree_path}()}
}
\concept{path finders}
