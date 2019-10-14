#' Intestinal microbiota of 84 individuals
#'
#' Species abundances of intestinal microbiota of 84 individuals, 41 with
#' Chron's Disease, 24 with Ulcerative Colitis, and 19 without inflamatory
#' inflammatory bowel disease.
#'
#' @format A data frame with 84 rows and 328 variables:
#' @source \url{https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products}
"speMGX"

#' Intestinal microbiota of 84 individuals
#'
#' Genus abundances of intestinal microbiota of 84 individuals, 41 with
#' Chron's Disease, 24 with Ulcerative Colitis, and 19 without inflamatory
#' inflammatory bowel disease.
#'
#' @format A data frame with 84 rows and 199 variables:
#' @source \url{https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products}
"genMGX"

#' Metagenome sampling times of 131 individuals (in weeks)
#'
#' @format A data frame with 131 rows (patients) and 55 columns (week number)
#' @source \url{https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products}
"visitTable"

#' Fecal Calprotectin levels
#'
#' @format A data frame with 131 rows (patients) and 52 columns (week number)
#' @source \url{https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv}
"FcalWeek"

#' Microbiota of 107 individuals
#'
#' @format A data frame with 1607 rows (patients) and 585 columns (species)
#' @source \url{https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products}
"allMGX"


#' Longitudinal data (intestinal microbiota)
#'
#' @format A list containing: a matrix with the species abundances of 107 individuals with 4 visits each one,
#' a vector with the number of individuals and visits, and a matrix with the distance among the visits (in weeks)
#' @source \url{https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products}
"longMGX"

#' Longitudinal metadata (intestinal microbiota)
#'
#' @format A table with the sampling time, individual diagnosis, sex, age etc
#' @source \url{https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products}
"metaLongMGX"
