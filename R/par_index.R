#' Protected area representativeness
#'
#' Using GDM predictions for assessing how well the protected areas in a
#' region represent the biological diversity.
#'
#' @param x A SpatRaster of GDM transformed layers. It accepts any raster type convertible to terra objects
#' or a path to files on disk.
#' @param pa protected areas as polygons, masked GDM transformed rasters, or a data.frame of extracted values
#'  from protected areas.
#' @param n number of samples to take from protected areas to compare with regions. NULL or 0 means all
#'  cells.
#' @param mod either a gdm model or the intercept of fitted model.
#' @param ncores number of cores to process. Default is all available cores.
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.
#'
#' @return SpatRaster
#' @export
#'
#' @examples
par_index <- function(x, pa, mod, n = NULL, ncores = -1, ...) {

    if (.is_rast(x)) {
        x <- .check_rast(x)
    } else {
        stop("x must be a raster object or path to raster objects.")
    }

    if (.is_vect(pa)) {
        pa_vals <- terra::extract(x, pa, ID = FALSE)
    } else if (.is_rast(pa)) {
        pa <- .check_rast(pa)
        pa_vals <- terra::as.data.frame(pa)
    } else if (methods::is(pa, "data.frame")) {
        pa_vals <- pa[, names(x)]
    } else {
        stop("pa must be a vector, raster or a data frame of GDM transform values for protected areas.")
    }

    n_ref <- nrow(pa_vals)
    if (is.null(n) || n == 0) {
        samps <- seq_len(n_ref) - 1 # to match C++ indices
    } else {
        samps <- seq_len(n_ref) - 1 # to match C++ indices
        if (n < n_ref) {
            samps <- sample(samps, n)
        }
    }

    out <- terra::predict(
        object = x,
        model = list(),
        fun = function(obj, dat, ...) {
            par_cpp(
                rast_vals = as.matrix(dat),
                ...
            )
        },
        ref_vals = as.matrix(pa_vals),
        samples = samps,
        intercept = ifelse(methods::is(mod, "gdm"), mod$intercept, mod),
        nthreads = ncores,
        ...
    )

    return(out)
}

