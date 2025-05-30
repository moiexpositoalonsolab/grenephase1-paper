
# GLM functions
glmfunctionp<-function(testdelta,testyear,testbio){
  lmmod<-lm(testdelta ~ testyear + testbio)
  p=coefficients(summary(lmmod))[3,4]
  b=coefficients(summary(lmmod))[3,1]
  return(p)
}
glmfunctionb<-function(testdelta,testyear,testbio){
  lmmod<-lm(testdelta ~ testyear + testbio)
  p=coefficients(summary(lmmod))[3,4]
  b=coefficients(summary(lmmod))[3,1]
  return(b)
}


# GLM binomial functions
glmfunctionpbinomial<-function(testfreq, testcount,testyear,testbio){
  lmmod<-glm(cbind(testfreq,testcount) ~ testyear + testbio,
             family = quasibinomial(link = logit))
  p=drop1(lmmod, test="Chisq")
  p=as.data.frame(p)[3,4]
  return(p)
}
