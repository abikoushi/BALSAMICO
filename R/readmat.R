readmat <- function(file, ...){
  dat <- scan(file,
              nlines=1, ...)
  ncol <- length(dat)
  dat <- matrix(
    scan(file),
    ncol = ncol, byrow = TRUE, ...)
  return(dat)
}
