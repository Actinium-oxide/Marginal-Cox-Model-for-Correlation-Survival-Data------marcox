\name{diabetes}
\docType{data}
\alias{diabetes}
\title{Diabetes Study Data}
\description{
  A dataset containing clinical information from a diabetes study.
}
\usage{data(diabetes)}
\format{
  A data frame with 166 rows and 6 variables:
  \describe{
    \item{\code{risk}}{Numeric: Risk score of the patient.}
    \item{\code{cens}}{Binary (0/1): Censoring indicator (1 = event occurred, 0 = censored).}
    \item{\code{time}}{Numeric: Time to event or censoring (in months).}
    \item{\code{id}}{Integer: Patient ID.}
    \item{\code{trt}}{Binary (0/1): Treatment indicator (1 = treated, 0 = control).}
    \item{\code{age}}{Binary (0/1): Age group indicator (1 = older, 0 = younger).}
  }
}
\source{
  Hypothetical clinical study data.
}
\examples{
  data(diabetes)
  summary(diabetes)
}
