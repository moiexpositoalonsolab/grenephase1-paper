# Title: Run linear mixed models in parallel
# Author: Meixi Lin
# Date: Thu Sep  7 12:02:10 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lmm_meixi")

library(dplyr)
library(data.table)
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

date()
sessionInfo()

# def functions --------
# assemble lmm objects
prep_lmm <- function(yy, env_sites, envvar, pop_strc) {
    mydata = cbind(yy, env_sites[,c('mergeid', 'site', envvar)], pop_strc)
    return(mydata)
}

# get lmm model results
format_lmm <- function(mymodel, envvar) {
    lmesum = summary(mymodel)
    lmer2 = MuMIn::r.squaredGLMM(mymodel) # contains two rows
    outdt = c(lmer2,
              lmesum$tTable[envvar, "Value"],
              lmesum$tTable[envvar, "p-value"],
              lmesum$BIC)
    return(outdt)
}

# def variables --------
# not batching, run full genome in one run
args = commandArgs(trailingOnly=TRUE)
#envvar = as.character(args[1])
#mygen = as.character(args[2])


envvar = 'bio1'
mygen = '1'
print(envvar)
print(mygen)

myfm = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:3), envvar), collapse = ' + ')))
print(myfm)

outdir = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lmm_meixi/results'
dir.create(outdir, recursive = TRUE)
dir.create(paste0(outdir, 'qqplots'))

# load data --------
load(paste0('./data/lmm_pc/inputs/envsites_popstrc_deltap_gen', mygen, '.rda'))
# deltap = deltap[1:1000, ]
print(dim(deltap))
# plot(pop_strc[,1:3])

# main --------
# run the linear mixed model in a parallel for loop
# note that this won't report the errors if there were any
# which means the result.ii row name is very meaningful if any error occurs 
lmeres <- foreach(ii = 1:nrow(deltap), .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    yy = as.numeric(unlist(deltap[ii,]))
    .GlobalEnv$myfm <- myfm # fix a global env bug
    mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
    model = nlme::lme(fixed = myfm, random = ~ 1|site, data = mydata) # 3 popstr PCs
    format_lmm(model, envvar) # output model results
}

dimnames(lmeres)[[2]] = c('R2m', 'R2c', 'beta', 'beta_p', 'BIC')

# these two number of rows can be different if there were errors during the lme model 
print(dim(deltap))
print(dim(lmeres))

# preliminary plotting
png(filename = paste0(outdir, 'qqplots/QQlme_', envvar, '_gen', mygen, '.png'), width = 10, height = 5, res = 300, units = 'in')
par(mfrow = c(1,2))
qqman::qq(lmeres[,'beta_p'], main = paste0('yy ~ PC1:3 + ', envvar, ' gen', mygen))
qqman::qq(lmeres[,'beta_p'], main = paste0('yy ~ PC1:3 + ', envvar, ' gen', mygen), xlim = c(0,7), ylim = c(0,7))
dev.off()

# output files --------
save(lmeres, file = paste0(outdir, 'lmeres_PC1to3_', envvar, '_gen', mygen, '.rda'))
# no RData image saved as everything is stored in lmeres or printed in the logs

# cleanup --------
stopCluster(cl)
date()
closeAllConnections()