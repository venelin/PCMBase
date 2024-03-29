% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{PCMPlotTraitData2D}
\alias{PCMPlotTraitData2D}
\title{Scatter plot of 2-dimensional data}
\usage{
PCMPlotTraitData2D(
  X,
  tree,
  sizePoints = 2,
  alphaPoints = 1,
  labeledTips = NULL,
  sizeLabels = 8,
  nudgeLabels = c(0, 0),
  palette = PCMColorPalette(PCMNumRegimes(tree), PCMRegimes(tree)),
  scaleSizeWithTime = !is.ultrametric(tree),
  numTimeFacets = if (is.ultrametric(tree) || scaleSizeWithTime) 1L else 3L,
  nrowTimeFacets = 1L,
  ncolTimeFacets = numTimeFacets
)
}
\arguments{
\item{X}{a k x N matrix}

\item{tree}{a phylo object}

\item{sizePoints, alphaPoints}{numeric parameters passed as arguments size and alpha to \code{geom_point}.
Default: sizePoints = 2, alphaPoints = 1.}

\item{labeledTips}{a vector of tip-numbers to label (NULL by default)}

\item{sizeLabels}{passed to \code{geom_text} to specify the size of tip-labels for the trait-points.}

\item{nudgeLabels}{a numeric vector of two elements (default: c(0,0)), passed as
arguments nudge_x and nudge_y of \code{geom_text}.}

\item{palette}{a named vector of colors}

\item{scaleSizeWithTime}{logical indicating if the size and the transparency of the points
should reflect the distance from the present (points that are farther away in time with
respect to the present moment, i.e. closer to the root of the tree, are displayed smaller
and more transparent.). By default this is set to \code{!is.ultrametric(tree)}.}

\item{numTimeFacets}{a number or a numeric vector controlling the creation of different facets
corresponding to different time intervals when the tree is non-ultrametric. If a single number,
it will be interpreted as an integer specifying the number of facets, each facets corresponding to
an equal interval of time. If a numeric vector, it will be used to specify the cut-points for
each interval. Default: \code{if(is.ultrametric(tree) || scaleSizeWithTime) 1L else 3}.}

\item{nrowTimeFacets, ncolTimeFacets}{integers specifying how the time facets should
be layed out. Default: \code{nrowTimeFacets = 1L, ncolTimeFacets = numTimeFacets}.}
}
\value{
a ggplot object
}
\description{
Scatter plot of 2-dimensional data
}
