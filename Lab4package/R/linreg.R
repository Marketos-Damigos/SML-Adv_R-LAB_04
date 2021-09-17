#' linreg RC Class
#'
#' @description This class offers linear regression calculation using the QR Decomposition method and comes with some handy tools.
#'
#' @field X matrix. Independent Values.
#' @field y matrix. Depenent Values.
#' @field b_hat matrix. Regression Coeffcients.
#' @field y_hat matrix. Fitted Values.
#' @field l_e matrix. Residuals.
#' @field l_n numeric. Number of Independant Values (Number of Observations).
#' @field l_p numeric. Number of Dependant Values (Number of Parameters in the Model).
#' @field l_df numeric. Degrees of freedom.
#' @field l_sigma_s matrix. Residual Variance.
#' @field l_var_beta matrix. Variance of the Regression Coefficients.
#' @field l_formula formula. The Formula for the Linear Regression.
#' @field l_p_values numeric. P-Values.
#' @field l_t_beta matrix. T-Values for each Coefficient.
#' @field l_data_set_name character. The Given Data.
#'
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
                      
                      # Methods ----------------------------
                      methods = list(
                        
                        # Constructor ----------------------
                        initialize = function(formula, data) {
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