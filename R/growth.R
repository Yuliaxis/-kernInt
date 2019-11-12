#' Growth data
#'
#' Longitudinal height data of 54 girls and 39 boys (93 individuals in total) from age 11 to 18.
#'
#' @format A matrix with 744 rows (93 individuals * 8 measures) and 3 variables:
#'  \describe{
#'   \item{Sex}{Boy or girl}
#'   \item{Age}{Ranging from 11 to 18}
#'   \item{Height}{Height in cm}
#' }
#' @source Berkeley Growth Study Data
#' Ramsay, James O., and Silverman, Bernard W. (2006), Functional Data Analysis, 2nd ed., Springer, New York.
#'
#' Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer, New York, ch. 6.
#'
#' Tuddenham, R. D., and Snyder, M. M. (1954) "Physical growth of California boys and girls from birth to age 18", University of California Publications in Child Development, 1, 183-364.
"growth"

#' Kernel matrices from Growth data
#'
#' @format A list with example kernel matrices from 'growth' data:
#'  \describe{
#'   \item{[[1]]}{Categorical kernel over the 'sex' column in growth.}
#'   \item{[[2]]}{RBF kernel over the 'age' and 'height' column (gamma=0.01)}
#'   \item{[[3]]}{RBF kernel over the 'age' column (gamma=0.05)}
#'   \item{[[4]]}{Number of individuals and number of measures}
#'   \item{[[5]]}{A table Individual x Age marking the available data with a 1 and the absent with a 0}
#' }
#' @source Berkeley Growth Study Data
#' Ramsay, James O., and Silverman, Bernard W. (2006), Functional Data Analysis, 2nd ed., Springer, New York.
#'
#' Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer, New York, ch. 6.
#'
#' Tuddenham, R. D., and Snyder, M. M. (1954) "Physical growth of California boys and girls from birth to age 18", University of California Publications in Child Development, 1, 183-364.
"grKern"
