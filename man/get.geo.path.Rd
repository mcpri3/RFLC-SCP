% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.geo.path.R
\name{get.geo.path}
\alias{get.geo.path}
\title{Border-to-border path within ecological continuities}
\usage{
get.geo.path(df, maxdisp, all.patches = all.ep)
}
\arguments{
\item{df}{a data.frame with a column from indicating the origin patch name and a column to indicating the destination patch(es) name.}

\item{maxdisp}{maximal path length (longer paths are removed).}

\item{all.patches}{sf object (polygons) of patch distribution. Patch name in all.patches should be indicated in a SITECODE column and match patch name in df.}
}
\value{
A data.frame with from-to (same as df) columns and a LinkLgth.km2 column indicating path length for paths that passed the filtering (i.e., <= maxdisp).
}
\description{
get.geo.path calculates the border-to-border shortest path between origin and destination patches belonging to the same ecological continuity.
}
