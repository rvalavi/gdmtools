# drop NA from file
drop_na <- function(x, war = TRUE, msg = ".") {
    n <- nrow(x)
    ccs <- stats::complete.cases(x)
    len <- sum(ccs)
    if (war && len < n) {
        message(n - len, " of records had NA and were removed", msg)
    }

    return(
        x[ccs, ]
    )
}

.is_mat <- function(x) {
    z <- class(x)
    return(
        any(
            z %in% c(
                "data.frame", "matrix", "data.table", "tibble"
            )
        )
    )
}


# is it a raster object
.is_rast <- function(x){
    z <- class(x)
    return(
        any(
            z %in% c(
                "SpatRaster",
                "RasterStack", "RasterLayer", "RasterBrick",
                "stars",
                "character"
            )
        )
    )
}


.is_vect <- function(x){
    z <- class(x)
    return(
        any(
            z %in% c(
                "SpatVector ", "sf"
            )
        )
    )
}

# is it a raster or convertible to raster?
.check_rast <- function(r, name = "x"){
    if(!methods::is(r, "SpatRaster")){
        tryCatch(
            {
                r <- terra::rast(r)
            },
            error = function(cond) {
                message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
                message(sprintf("'%s' must be a SpatRaster, stars, Raster* object, or path to a raster file on disk.", name))
            }
        )
    }
    return(r)
}

# check for required packages
.check_pkgs <- function(pkg){
    pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
    if(length(pkgna) > 0){
        nm <- paste(pkgna, collapse = ", ")
        message("This function requires these packages: ", nm, "\nWould you like to install them now?\n1: yes\n2: no")
        user <- readline(prompt = paste0("Selection: "))
        if(tolower(user) %in% c("1", "yes", "y")){
            utils::install.packages(pkgna)
        } else{
            stop("Please install these packages for function to work: ", nm)
        }
    }
}
