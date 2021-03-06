% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treatmentResponseDiallel.R
\name{var.translate}
\alias{var.translate}
\title{var.translate: translate variable names}
\usage{
var.translate(variable.name, ...)
}
\arguments{
\item{variable.name}{this is the name or names of variables to be translated}

\item{...}{additional arguments}
}
\value{
returns translation result
}
\description{
Translate variable names to simpler/more readable standards.
}
\examples{
strings <- c("tau:Gender:aj", "SymCrossjk", ":j:4")
var.translate(strings)
}
