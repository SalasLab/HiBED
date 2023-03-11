#' HiBED library CpGs matrix for brain tissue DNA methylation deconvolution
#'
#' @description
#'     This object contains 4 matrices of the the average
#'     DNA methylation values of the probes included in 4 layers of the HiBED
#'     deconvolution. These CpGs are used as the backbone for deconvolution and
#'     were selected because their methylation signature differs across the
#'     seven brain cell subtypes.
#'
#' @format The matrices are 81 x 3, 183 x 4, 237 x 5, 120 x 4
#'
#'         The format is:
#'         num [1:81, 1:3] 0.04592944  0.02268472  0.88886150 ...
#'
#' @examples
#' # data ("HiBED_Libraries")
#' # head(Library_Layer1)
"HiBED_Libraries"
