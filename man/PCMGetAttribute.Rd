% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCM.R
\name{PCMGetAttribute}
\alias{PCMGetAttribute}
\title{Value of an attribute of an object or values for an attribute found in its members}
\usage{
PCMGetAttribute(name, object, member = "", ...)
}
\arguments{
\item{name}{attribute name.}

\item{object}{a PCM model object or a PCMTree object.}

\item{member}{a member expression. Member expressions are character strings
denoting named elements in a list object (see examples). Default: "".}

\item{...}{additional arguments passed to \code{\link{MatchListMembers}}.}
}
\value{
if member is an empty string, \code{attr(object, name)}. Otherwise, a named list
containing the value for the attribute for each member in \code{object}
matched by \code{member}.
}
\description{
Value of an attribute of an object or values for an attribute found in its members
}
\examples{
PCMGetAttribute("class", PCMBaseTestObjects$model_MixedGaussian_ab)
PCMGetAttribute(
  "dim", PCMBaseTestObjects$model_MixedGaussian_ab,
  member = "$Sigmae_x")
}
