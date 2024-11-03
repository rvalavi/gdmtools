#' Protected area representativeness
#'
#' Using GDM predictions for assessing how well the protected areas in a
#' region represent the biological diversity.
#'
#' Make sure \href{https://en.wikipedia.org/wiki/OpenMP}{OpenMP} is installed on your system to enable
#' parallel processing and speed up computations.
#'
#' @param x A SpatRaster of GDM-transformed layers. Accepts any raster type that can be converted to a \pkg{terra} object
#' or a file path to rasters stored on disk.
#' @param pa Protected areas provided as polygons, masked GDM-transformed rasters, or a data.frame of values
#' extracted GDM transform grids from protected areas.
#' @param n A numeric value specifying the number of samples to take from the entire region
#' to normalize similarity values. If set to \code{NULL} or 0, 10\% of the non-NA cells will be used
#' as the sample size. Alternatively, a data.frame of extracted values from random cells within the
#' region can be provided instead.
#' @param mod Either a \pkg{gdm} model object or the intercept value from a fitted model.
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
par_index <- function(x, pa, mod, n = NULL, ncores = -1, ...) {

    if (.is_rast(x)) {
        x <- .check_rast(x)
    } else {
        stop("x must be a raster object or path to raster objects.")
    }

    # get the pa values
    pa_vals <- {
        if (.is_vect(pa)) {
            as.matrix(terra::extract(x, pa, ID = FALSE))
        } else if (.is_rast(pa)) {
            as.matrix(terra::as.data.frame(.check_rast(pa)))
        } else if (.is_mat(pa)) {
            as.matrix(pa[, names(x)])
        } else {
            stop("pa must be a vector, raster or a data frame of GDM transform values for protected areas.")
        }
    };

    # filter NAs
    pa_vals <- drop_na(pa_vals, msg = " from the protected area samples.")

    rand_cells <- {
        if (.is_mat(n)) {
            as.matrix(n)
        } else {
            # take random sample from the region
            if (is.null(n) || n < 1) {
                nc <- as.numeric(terra::global(x[[1]], fun = "notNA"))
                n <- ceiling(nc / 10)
            }
            as.matrix(
                drop_na(
                    terra::spatSample(x, size = n, method = "random"),
                    msg = " from the random samples."
                )
            )
        }
    };

    out <- terra::predict(
        object = x,
        model = list(),
        fun = function(model, newdata, ...) {
            return(
                par_cpp(
                    rast_vals = as.matrix(newdata),
                    ...
                )
            )
        },
        ref_vals = pa_vals,
        samp_vals = rand_cells,
        intercept = ifelse(methods::is(mod, "gdm"), mod$intercept, mod),
        nthreads = ncores,
        ...
    )

    return(out)
}

# # wrapper function for par C++
# par_fun <- function(model, newdata, ...) {
#     return(
#         par_cpp(
#             rast_vals = as.matrix(newdata),
#             ...
#         )
#     )
# }

