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

## load metadata
env_sites = data.table::fread(file = 'metadt_first_gen.csv')
## load pop structure
pop_strc = data.table::fread(file = 'pop_structure_first_gen.csv')
#delta_p_filt.to_csv('../key_files/delta_p_maf05_firstgensamples.csv',index=None)
infile = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/key_files/delta_p_maf05_mincount05_firstgensamples.csv'
#key_files/delta_p_maf05_firstgensamples.csv'
## load delta p 
## the selec =-1 is because the frist columns is the idnex 
deltap = data.table::fread(file = infile, drop = 1)
print(dim(deltap))

# def functions --------
# assemble lmm objects
#prep_lmm <- function(yy, env_sites, pop_strc) {
#  mydata <- yy %>%
#    bind_cols(env_sites %>% select(mergeid, site, bio1), pop_strc)
#  return(mydata)
#}
prep_lmm <- function(yy, env_sites, envvar, pop_strc) {
    mydata = cbind(yy, env_sites[,c('mergeid', 'site', 'bio1')], pop_strc)
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
#args = commandArgs(trailingOnly=TRUE)
#envvar = as.character(args[1])
#mygen = as.character(args[2])
mygen = 'first_gen'
envvar = 'bio1'
print(envvar)

myfm = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:3), envvar), collapse = ' + ')))
print(myfm)

outdir = '/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lmm_meixi/results/'
dir.create(outdir, recursive = TRUE)
dir.create(paste0(outdir, 'qqplots'))


# main --------
# run the linear mixed model in a parallel for loop
# note that this won't report the errors if there were any
# which means the result.ii row name is very meaningful if any error occurs 
lmeres <- foreach(ii = 1:nrow(deltap), .combine = 'rbind', .errorhandling = 'pass') %dopar% {
  tryCatch({
    yy = as.numeric(unlist(deltap[ii,]))
    .GlobalEnv$myfm <- myfm # fix a global env bug
    mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
    model = nlme::lme(fixed = myfm, random = ~ 1|site, data = mydata) # 3 popstr PCs
    format_lmm(model, envvar) # output model results
  }, error = function(e) {
    # Log or handle the error
    message(paste("Error at row", ii, ":", e$message))
    return(NA) # return a value (e.g., NA) if an error occurs
  })
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
#save(lmeres, file = paste0(outdir, 'lmeres_PC1to3_', envvar, '_gen', mygen, '.rda'))
write.csv(lmeres, file = paste0(outdir, 'lmeres_PC1to3_', envvar, '_gen', mygen, '.csv'), row.names = FALSE)

# no RData image saved as everything is stored in lmeres or printed in the logs

# cleanup --------
stopCluster(cl)
date()
closeAllConnections()