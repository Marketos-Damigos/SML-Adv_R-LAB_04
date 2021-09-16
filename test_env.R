linreg <- function(formula, data){
  X = model.matrix(formula, data)
  y = as.matrix(data[all.vars(formula)[1]])
  QR = qr(X)
  Q = qr.Q(QR)
  R = qr.R(QR)
  Qy = t(Q) %*% y
  b_hat = solve(R) %*% Qy
  y_hat = round(Q %*% t(Q) %*% y, 5)
  e_hat = round(y - y_hat, 5)
  df = length(X[,1]) - length(X[1,])
  s2_hat = (t(e_hat) %*% e_hat) / df
  var_b_hat = as.numeric(s2_hat) * solve(crossprod(R))
  t_b = b_hat / sqrt(diag(var_b_hat))
  p_value = 2*pt(-abs(t_b), df)
  
  residuals_vs_fitted <- ggplot(data.frame(x = e_hat, y = y_hat), aes(y_hat, e_hat)) +
    geom_point() +
    xlab(paste("Fitted Values\n", "linreg(", format(formula), ")", ""))+
    stat_summary(fun = median, colour="#54d8e0", geom="line",group=1)
  
  scale_location = ggplot(stdresfit, aes(x = fitted, y = stdResiduals))+
    geom_point(color = "white")+
    xlab(paste("Fitted Values\n", "lm(", format(l_formula), ")", ""))+
    ylab(expression(sqrt("|Standardized residuals|")))+
    stat_summary(aes(y = stdResiduals, x = fitted ,group=1),
                 fun= mean, colour="#54d8e0", geom="line",group=1) + liu_theme()

  library(plotly)
  z = ggplotly(residuals_vs_fitted)

}

linreg(formula = Petal.Length~Species, data = iris)


