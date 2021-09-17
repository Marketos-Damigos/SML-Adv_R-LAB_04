asterisks <- function(p_value) {
    if (p_value > 0.1) return(" ")
    if (p_value > 0.05) return(".")
    if (p_value > 0.01) return("*")
    if (p_value > 0.001) return("**")
    return("***")
}

asterisks <- Vectorize(asterisks, vectorize.args = "p_value")