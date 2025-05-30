lm_eq<-
function (y, x, tex = TRUE)
{
  mylm <- lm(y ~ x)
  p = format(coefficients(summary(mylm))[2, 4], digits = 3,
             scientific = TRUE)
  r2 = round(summary(mylm)$r.squared, digits = 3)
  b = format(coefficients(summary(mylm))[2, 1], digits = 3,
            scientific=TRUE)
  if (tex == FALSE) {
    return(sprintf("R2= %s, b= %s, p= %s", r2, b, p))
  }
  else {
    return(paste0("$R^2 = $", r2, ", $\\beta = $", b, ", $ p = $",
                  p))
  }
}

ht<-function(d){
  rbind(head(d), tail(d))
}
