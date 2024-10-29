#' Protected area representativeness
#'
#' Using GDM predictions for assessing how well the protected areas in a
#' region represent the biological diversity.
#'
#' @param x A SpatRaster of GDM transformed layers. It accepts any raster type convertible to terra objects
#' or a path to files on disk.
#' @param pa protected areas as polygons, masked GDM transformed rasters, or a data.frame of extracted values
#'  from protected areas.
#' @param n number of samples to take from the entire region to normalise the similarity values. if NULL or 0
#' used, 20% of the number of non-na cells will be used as the sample number.
#' @param mod either a gdm model or the intercept of fitted model.
#' @param ncores number of cores to process. Default is all available cores.
#' @param ... Additional arguments for writing raster outputs e.g. \code{filename},
#' \code{overwrite}, and \code{wopt} from terra \code{\link[terra]{predict}}.
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' library(gdmtools)
#'
par_index <- function(x, pa, mod, n = NULL, ncores = -1, ...) {

    if (.is_rast(x)) {
        x <- .check_rast(x)
    } else {
        stop("x must be a raster object or path to raster objects.")
    }

    if (.is_vect(pa)) {
        pa_vals <- as.matrix(terra::extract(x, pa, ID = FALSE))
    } else if (.is_rast(pa)) {
        pa <- .check_rast(pa)
        pa_vals <- as.matrix(terra::as.data.frame(pa))
    } else if (methods::is(pa, "data.frame")) {
        pa_vals <- as.matrix(pa[, names(x)])
    } else {
        stop("pa must be a vector, raster or a data frame of GDM transform values for protected areas.")
    }

    # filter NAs
    pa_vals <- pa_vals[stats::complete.cases(pa_vals), ]

    if (is.null(n) || n == 0) {
        nc <- as.numeric(terra::global(x[[1]], fun = "notNA"))
        n <- ceiling(nc / 5)
    }
    rand_cells <- terra::spatSample(x, size = n, method = "random")
    rand_cells <- as.matrix(rand_cells[stats::complete.cases(rand_cells), ])

    par_fun <- function(model, newdata, ...) {
        # check for NAs
        has_na <- anyNA(newdata)
        nr <- nrow(newdata)
        if (has_na) {
            idx <- which(stats::complete.cases(newdata))
            out <- rep(NaN, nr)
            # if all NA, return NaN vector
            if (!length(idx)) return(out)
            # subset the complete data
            dat <- as.matrix(newdata[idx, ])
        } else {
            dat <- as.matrix(newdata)
        }

        pcount <- par_cpp(
            rast_vals = dat,
            ...
        )

        # sort out possible NAs
        if (has_na) {
            out[idx] <- pcount
            return(out)
        } else {
            return(pcount)
        }
    }

    out <- terra::predict(
        object = x,
        model = list(),
        fun = par_fun,
        ref_vals = pa_vals,
        samp_vals = rand_cells,
        intercept = ifelse(methods::is(mod, "gdm"), mod$intercept, mod),
        nthreads = ncores,
        ...
    )

    return(out)
}

