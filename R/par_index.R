#' Protected area representativeness
#'
#' Using GDM predictions for assessing how well the protected areas in a
#' region represent the biological diversity.
#'
#' @param x
#' @param pa
#' @param n
#' @param intercept
#' @param ncores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
par_index <- function(x, pa, intercept, n = NULL, ncores = -1, ...) {

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
        samps <- seq_len(n_ref)
    } else {
        n <- ifelse(n > n_ref, n_ref, n)
        samps <- seq_len(n_ref)
    }

    out <- terra::predict(
        object = r,
        model = list(),
        fun = function(obj, dat, ...) {
            par_cpp(
                rast_vals = as.matrix(dat),
                ...
            )
        },
        ref_vals = as.matrix(pa_vals),
        samples = samps,
        intercept = ifelse(methods::is(intercept, "gdm"), intercept$intercept, intercept),
        nthreads = ncores,
        ...
    )

    return(out)
}

