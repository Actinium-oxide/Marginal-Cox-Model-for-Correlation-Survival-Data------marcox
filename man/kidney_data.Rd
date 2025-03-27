\name{kidney_data}
\docType{data}
\alias{kidney_data}
\title{Kidney Disease Study Data}
\description{
  A dataset containing survival analysis information related to kidney disease patients.
}
\usage{data(kidney_data)}
\format{
  A data frame with 76 rows and 5 variables:
  \describe{
    \item{\code{time}}{Numeric: Time to event or censoring (in days).}
    \item{\code{cens}}{Binary (0/1): Censoring indicator (1 = event occurred, 0 = censored).}
    \item{\code{age}}{Numeric: Age of the patient in years.}
    \item{\code{sex}}{Binary (0/1): Sex of the patient (1 = male, 0 = female).}
    \item{\code{type}}{Categorical (0,1,2,3): Kidney disease type classification.}
  }
}
\source{
  Hypothetical survival study data.
}
\examples{
  data(kidney_data)
  summary(kidney_data)
}
