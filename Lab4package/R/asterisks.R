#' Asterisks function
#'
#' @description Help function to convert p_value to asterisks
#' @param p_value p value
#' @usage asterisks(p_value)
#' @return Returns the p_value in asterisks
#' @export asterisks



asterisks <- function(p_value) {
    if (p_value > 0.1) return(" ")
    if (p_value > 0.05) return(".")
    if (p_value > 0.01) return("*")
    if (p_value > 0.001) return("**")
    return("***")
}

asterisks <- Vectorize(asterisks, vectorize.args = "p_value")