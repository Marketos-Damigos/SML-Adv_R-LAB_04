#' linreg RC Class
#'
#' @description This class offers linear regression calculation using the QR Decomposition method and comes with some handy tools.
#'
#' @field X matrix. Independent values.
#' @field y matrix. Dependent values.
#' @field b_hat matrix. Regression coefficients.
#' @field y_hat matrix. Fitted values.
#' @field e_hat matrix. Residuals values.
#' @field df numeric. Degrees of freedom.
#' @field s2_hat matrix. Variance of residuals.
#' @field var_b_hat matrix. Variance of the regression coefficients.
#' @field f_text formula. The formula for the linear regression.
#' @field p_value numeric. P-Values.
#' @field t_b matrix. T-Values for each coefficient.
#' @field name_ds character. The given data.
#' @field Q matrix. Q matrix from the QR decomposition.
#' @field R matrix. Upper triangular matrix.
#' @field Qy matrix. 
#' @field sterr matrix. Standard error.
#' @field print() function. Returns coefficients.
#' @field plot function. Returns sets of plots.. 
#' @field resid function. Returns residuals. 
#' @field pred function. Returns predictions.
#' @field coef function. Returns coefficients.
#' @field summary function. Returns summary of the above.
#' @import methods
#' @return Nothing.
#' @export linreg
#' @exportClass linreg


linreg <- setRefClass("linreg",
                      fields = list(
                        X = "matrix",
                        y = "matrix",
                        b_hat = "matrix",
                        y_hat = "matrix",
                        e_hat = "matrix",
                        df = "numeric",
                        s2_hat = "matrix",
                        var_b_hat = "matrix",
                        sterr = "matrix",
                        f_text = "formula",
                        p_value = "matrix",
                        t_b = "matrix",
                        name_ds = "character",
                        Q = "matrix",
                        R = "matrix",
                        Qy = "matrix"),
                      
                      methods = list(
                        
                        initialize = function(formula, data) {
                          
                          if (class(formula) != "formula") {
                            stop("Error: You should give a formula")
                          }
                          if (is.data.frame(data) == FALSE){
                            stop("Error: You should give a dataframe")
                          }
                          
                          X <<- model.matrix(formula, data)
                          y <<- as.matrix(data[all.vars(formula)[1]])
                          QR <- qr(X)
                          Q <<- qr.Q(QR)
                          R <<- qr.R(QR)
                          Qy <<- t(Q) %*% y
                          b_hat <<- solve(R) %*% Qy
                          y_hat <<- round(Q %*% t(Q) %*% y, 5)
                          e_hat <<- round(y - y_hat, 5)
                          df <<- length(X[,1]) - length(X[1,])
                          s2_hat <<- (t(e_hat) %*% e_hat) / df
                          var_b_hat <<- as.numeric(s2_hat) * solve(crossprod(R))
                          sterr <<- abs((e_hat - mean(e_hat)))
                          t_b <<- b_hat / sqrt(diag(var_b_hat))
                          p_value <<- 2*pt(-abs(t_b), df)
                          
                          f_text <<- formula
                          name_ds <<- deparse(substitute(data))
                        },
                        

                        print = function(){
                          cat(paste("linreg(formula = ", format(f_text), ", data = ", name_ds, ")\n\n", sep = ""))
                          cat(paste("Coefficients:\n\n"))

                          ext_print(t(b_hat))
                        },


                        plot = function(){
                          
                          library(ggplot2)
                          
                          residuals_vs_fitted <- ggplot(data.frame(x = e_hat, y = y_hat), aes(y_hat, e_hat)) +
                            geom_point() +
                            xlab(paste("Fitted Values\n", "linreg(", format(formula), ")", ""))+
                            stat_summary(fun = median, colour="#54d8e0", geom="line",group=1)

                          ext_print(residuals_vs_fitted)

                          scale_location = ggplot(data.frame(x = sterr, y = y_hat), aes(x = y_hat, y = sterr))+
                            geom_point()+
                            xlab(paste("Fitted Values\n", "lm(", format(formula), ")", ""))+
                            ylab(expression(sqrt("|Standardized residuals|")))+
                            stat_summary(fun= mean, colour="#54d8e0", geom="line",group=1)

                          ext_print(scale_location)
                        },

                        resid = function() {
                          return(e_hat)
                        },

                        pred = function() {
                          return(y_hat)
                        },

                        coef = function() {
                          return(b_hat)
                        },

                        summary = function() {
                          
                          cat("\nCall:\n")
                          cat(paste("linreg(formula = ", (format(f_text)), ", data = ", name_ds, ")\n\n", sep = ""))
                          

                          cat("Coefficients:\n")
                          
                          sum_table = data.frame(b_hat, sqrt(diag(var_b_hat)), t_b, p_value, asterisks(p_value))
                          colnames(sum_table) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
                          ext_print(sum_table)
                          cat(paste("\nResidual standard error:", sqrt(s2_hat[1]), "on", df, "degrees of freedom"))
                        }
                        
                      )
)