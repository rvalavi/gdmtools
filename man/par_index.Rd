% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/par_index.R
\name{par_index}
\alias{par_index}
\title{Protected area representativeness}
\usage{
par_index(x, pa, mod, n = NULL, ncores = -1, ...)
}
\arguments{
\item{x}{A SpatRaster of GDM-transformed layers. Accepts any raster type that can be converted to a \pkg{terra} object
or a file path to rasters stored on disk.}

\item{pa}{Protected areas provided as polygons, masked GDM-transformed rasters, or a data.frame of values
extracted GDM transform grids from protected areas.}

\item{mod}{Either a \pkg{gdm} model object or the intercept value from a fitted model.}

\item{n}{A numeric value specifying the number of samples to take from the entire region
to normalize similarity values. If set to \code{NULL} or 0, 10\% of the non-NA cells will be used
as the sample size. Alternatively, a data.frame of extracted values from random cells within the
region can be provided instead.}

\item{ncores}{Integer. Specifies the number of CPU threads to be used for processing. A value
below 1 indicates that all available threads will be utilized. Refer to the details section for
more information.}

\item{...}{Additional arguments for writing raster outputs e.g. \code{filename},
\code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.}
}
\value{
SpatRaster
}
\description{
Using GDM predictions for assessing how well the protected areas in a
region represent the biological diversity.
}
\details{
Make sure \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system to enable
parallel processing and speed up computations.
}
\examples{
# library(gdm)
library(gdmtools)

}
\references{
Mokany, K., Ware, C., Woolley, S. N. C., Ferrier, S., & Fitzpatrick, M. C. (2022).
A working guide to harnessing generalized dissimilarity modelling for biodiversity analysis and
conservation assessment. Global Ecology and Biogeography, 31(4), 802–821.
}
