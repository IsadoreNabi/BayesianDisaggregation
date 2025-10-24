\name{save_results}
\alias{save_results}
\title{Save disaggregation results to disk}
\usage{
save_results(res, out_dir)
}
\arguments{
\item{res}{A result object returned by \code{bayesian_disaggregate()}.}
\item{out_dir}{Directory where files will be written. Created if missing.}
}
\value{
(Invisibly) the path to the Excel file written.
}
\description{
Writes CSV extracts and a single Excel workbook with the key outputs from
\code{\link{bayesian_disaggregate}}.
}
\examples{
\donttest{
# Usar archivos de ejemplo incluidos en el paquete
cpi_file <- system.file("extdata", "CPI.xlsx", 
                         package = "BayesianDisaggregation")
weights_file <- system.file("extdata", "WEIGHTS.xlsx", 
                            package = "BayesianDisaggregation")

if (nzchar(cpi_file) && nzchar(weights_file)) {
  res <- bayesian_disaggregate(cpi_file, weights_file)
  # Use tempdir() for CRAN compliance
  save_results(res, file.path(tempdir(), "results"))
}
}
}