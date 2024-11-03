#' Estimate Expected Species Persistence based on Habitat Condition
#'
#' This function combines spatial predictions from Generalized Dissimilarity
#' Modelling (GDM) with habitat condition data to estimate the expected level
#' of species persistence at each location.
#'
#' @param x A SpatRaster of GDM-transformed layers. Accepts any raster type that can be converted
#' to a \pkg{terra} object or a file path to rasters stored on disk.
#' @param y A SpatRaster of habitat condition. Accepts any raster type that can be converted
#' to a \pkg{terra} object or a file path to rasters stored on disk.
#' @param mod Either a \pkg{gdm} model object or the intercept value from a fitted model.
#' @param power Numeric. The power of z for species-area conversion (default 0.25).
#' @param n Numeric. Specifies the number of samples to take from the GDM transformed layers
#' and habitat condition data. If \code{NULL} or 0, a sample size of 10\% of the non-NA cells
#' will be used. This random sample may include NA values, which will be removed. Consequently,
#' the final number of samples used for analysis may be lower than the specified sample size.
#' @param ncores Integer. Specifies the number of CPU threads to be used for processing. A value
#' below 1 indicates that all available threads will be utilized. Refer to the details section for
#' more information.
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.
#'
#' @references Mokany, K., Ware, C., Woolley, S. N. C., Ferrier, S., & Fitzpatrick, M. C. (2022).
#' A working guide to harnessing generalized dissimilarity modelling for biodiversity analysis and
#' conservation assessment. Global Ecology and Biogeography, 31(4), 802–821.
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' # library(gdm)
#' library(gdmtools)
#'
persistance <- function(x, y, mod, power = 0.25, n = NULL, ncores = -1, ...) {

    if (.is_rast(x)) {
        x <- .check_rast(x, "x")
    } else {
        stop("x must be a raster object or path to raster objects.")
    }

    if (.is_rast(y)) {
        y <- .check_rast(y, "y")
    } else {
        stop("y must be a raster object or path to raster objects.")
    }


    # take random sample from the region
    if (is.null(n) || n < 1) {
        nc <- as.numeric(terra::global(x[[1]], fun = "notNA"))
        n <- ceiling(nc / 10)
    }
    rand_cells <- as.matrix(
        drop_na(
            terra::spatSample(c(x, y[[1]]), size = n, method = "random"),
            msg = " from the random samples."
        )
    )

    out <- terra::predict(
        object = x,
        model = list(),
        fun = function(model, newdata, ...) {
            return(
                presist_cpp(
                    rast_vals = as.matrix(newdata),
                    ...
                )
            )
        },
        ref_vals = rand_cells[, -ncol(rand_cells)],
        cond_vals = rand_cells[, ncol(rand_cells), drop = FALSE],
        intercept = ifelse(methods::is(mod, "gdm"), mod$intercept, mod),
        power = power,
        nthreads = ncores,
        ...
    )

    return(out)
}
