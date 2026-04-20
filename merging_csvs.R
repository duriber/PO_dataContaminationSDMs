## script to merge all csv saved model evaluations and coefficients into a single csv.

list_results <- list.files(path = 'out',
                           pattern = 'results',
                           full.names = TRUE)

# identifying the master file: "results.csv"
m <- which(lapply(list_results, function(x)
  {length(unlist(strsplit(x, '_')))}) ==1)

results_master_ISDM <- read.csv(list_results[m])
results_master_PPM_noBias <- read.csv(list_results[m])
results_master_PPM_withBias <- read.csv(list_results[m])

# list results ISDMs
list_results_ISDM <- list_results[which(lapply(list_results, function(x)
{length(unlist(strsplit(x, '_')))}) ==3)] 

for(i in 1:length(list_results_ISDM)){
  file_i <- list_results_ISDM[i]
  results_i <- read.csv(file_i)
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- unlist(strsplit(temp, '_'))[2]
  iter <-  unlist(strsplit(temp, '_'))[3]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  columns <- grep(model, x = names(results_i))
  results_master_ISDM[iter, columns] <- results_i[iter, columns]
}

write.csv(results_master_ISDM, file = 'out/results_master_ISDM.csv',
          row.names = TRUE)

# list results PPMs no Bias
list_results_PPM_noBias <- list.files(path = 'out',
                                      pattern = 'results_ppm_noBias',
                                      full.names = TRUE)


for(i in 1:length(list_results_PPM_noBias)){
  file_i <- list_results_PPM_noBias[i]
  results_i <- read.csv(file_i)
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- unlist(strsplit(temp, '_'))[4]
  iter <-  unlist(strsplit(temp, '_'))[5]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  columns <- grep(model, x = names(results_i))
  results_master_PPM_noBias[iter, columns] <- results_i[iter, columns]
}

write.csv(results_master_PPM_noBias, file = 'out/results_master_PPM_noBias.csv',
          row.names = TRUE)

# list results PPMs WITH Bias-related covariates
list_results_PPM_withBias <- list.files(path = 'out',
                                      pattern = 'results_ppm_withBias',
                                      full.names = TRUE)

for(i in 1:length(list_results_PPM_withBias)){
  file_i <- list_results_PPM_withBias[i]
  results_i <- read.csv(file_i)
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- unlist(strsplit(temp, '_'))[4]
  iter <-  unlist(strsplit(temp, '_'))[5]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  columns <- grep(model, x = names(results_i))
  results_master_PPM_withBias[iter, columns] <- results_i[iter, columns]
}

write.csv(results_master_PPM_withBias, file = 'out/results_master_PPM_withBias.csv',
          row.names = TRUE)


################################################################################
# now for the coefficients for ISDMs
list_coefs <- list.files(path = 'out',
                           pattern = 'coefficients',
                           full.names = TRUE)

# identifying the master file: "coefs_master.csv"
n <- which(lapply(list_coefs, function(x)
{length(unlist(strsplit(x, '_')))}) ==1)

coefs_master_ISDM <- 
  coefs_master_PPM_noBias <-
  coefs_master_PPM_withBias <-
  read.table(list_coefs[n],
             header = TRUE, sep = ",", row.names = 1)


list_coefs_ISDM <- list.files(path = 'out',
                              pattern = 'coefficients_m',
                              full.names = TRUE)


for(i in 1:length(list_coefs_ISDM)){
  file_i <- list_coefs_ISDM[i]
  coefs_i <- read.table(file_i,
                        header = TRUE, sep = ",", row.names = 1)
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- paste0(unlist(strsplit(temp, '_'))[2], '_')
  iter <-  unlist(strsplit(temp, '_'))[3]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))+2
  rows <- grep(model, x = rownames(coefs_i))
  coefs_master_ISDM[rows, iter] <- coefs_i[rows, iter]
}

write.csv(coefs_master_ISDM, 'out/coefs_master_ISDM.csv',
          row.names = TRUE)




list_coefs_PPM_noBias <- list.files(path = 'out',
                              pattern = 'coefficients_ppm_noBias',
                              full.names = TRUE)

for(i in 1:length(list_coefs_PPM_noBias)){
  file_i <- list_coefs_PPM_noBias[i]
  coefs_i <- read.table(file_i,
                        header = TRUE, sep = ",", row.names = 1)
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- paste0(unlist(strsplit(temp, '_'))[4], '_')
  iter <-  unlist(strsplit(temp, '_'))[5]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))+2
  rows <- grep(model, x = rownames(coefs_i))
  coefs_master_PPM_noBias[rows, iter] <- coefs_i[rows, iter]
}

write.csv(coefs_master_PPM_noBias, 'out/coefs_master_PPM_noBias.csv',
          row.names = TRUE)



list_coefs_PPM_withBias <- list.files(path = 'out',
                              pattern = 'coefficients_ppm_withBias',
                              full.names = TRUE)

for(i in 1:length(list_coefs_PPM_withBias)){
  file_i <- list_coefs_PPM_withBias[i]
  coefs_i <- read.table(file_i,
                        header = TRUE, sep = ",", row.names = 1)
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- paste0(unlist(strsplit(temp, '_'))[4], '_')
  iter <-  unlist(strsplit(temp, '_'))[5]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))+2
  rows <- grep(model, x = rownames(coefs_i))
  coefs_master_PPM_withBias[rows, iter] <- coefs_i[rows, iter]
}

write.csv(coefs_master_PPM_withBias, 'out/coefs_master_PPM_withBias.csv',
          row.names = TRUE)


################################################################################
# Finally for model predictions: evalDF CSVs

# first we create the 30 fataframes (1 per replicate) based on the evaluation dataset
preds_RDFs <- list.files(path = 'out/', 
                         pattern = 'PA_eval', 
                         full.names = TRUE)
for(k in 1:100){
  file_k <- preds_RDFs[k]
  temp <- tail(unlist(strsplit(file_k, '/')), n = 1)
  iter <-  unlist(strsplit(temp, '_'))[3]
  iter <- as.numeric(unlist(strsplit(iter, '.rds')))
  SimDat_i <- readRDS(file_k)
  write.csv(as.data.frame(SimDat_i)[,c('x', 'y', 'PA')], 
            file = paste0('out/Eval_iter_', iter, '.csv'),
            row.names = FALSE)
}

preds_iter <- list.files(path = 'out/', 
                         pattern = 'Eval_iter_', 
                         full.names = TRUE)


# for ISDMs

preds_CSVs_ISDM <- list.files(path = 'out/', 
                         pattern = 'evalDF_m', 
                         full.names = TRUE)

for(i in 1:length(preds_CSVs_ISDM)){
  file_i <- preds_CSVs_ISDM[i]
  preds_i <- read.csv(file_i,
                        header = TRUE)[,5:9]
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- paste0(unlist(strsplit(temp, '_'))[2], '_')
  model <- paste0("ISDM_", model)
  names(preds_i) <- paste0(model, names(preds_i))
  iter <-  unlist(strsplit(temp, '_'))[3]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  DF_i <- read.csv(preds_iter[grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)])
  DF_i <- cbind(DF_i, preds_i[,c(setdiff(names(preds_i), names(DF_i)))])
  DF_i <- write.csv(DF_i, 
                    file = preds_iter[grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)],
                    row.names = FALSE)
}


# for PPMs with no bias-related covariates

preds_CSVs_PPM_noBias <- list.files(path = 'out/', 
                              pattern = 'evalDF_ppm_noBias', 
                              full.names = TRUE)

for(i in 1:length(preds_CSVs_PPM_noBias)){
  file_i <- preds_CSVs_PPM_noBias[i]
  preds_i <- read.csv(file_i,
                      header = TRUE)[,5:9]
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- paste0(unlist(strsplit(temp, '_'))[4], '_')
  model <- paste0("PPM_noBias_", model)
  names(preds_i) <- paste0(model, names(preds_i))
  iter <-  unlist(strsplit(temp, '_'))[5]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  DF_i <- read.csv(preds_iter[grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)])
  DF_i <- cbind(DF_i, preds_i[,c(setdiff(names(preds_i), names(DF_i)))])
  DF_i <- write.csv(DF_i, 
                    file = preds_iter[grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)],
                    row.names = FALSE)
}


# for PPMs with bias-related covariates

preds_CSVs_PPM_withBias <- list.files(path = 'out/', 
                                    pattern = 'evalDF_ppm_withBias', 
                                    full.names = TRUE)

for(i in 1:length(preds_CSVs_PPM_withBias)){
  file_i <- preds_CSVs_PPM_withBias[i]
  preds_i <- read.csv(file_i,
                      header = TRUE)[,5:9]
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- paste0(unlist(strsplit(temp, '_'))[4], '_')
  model <- paste0("PPM_withBias_", model)
  names(preds_i) <- paste0(model, names(preds_i))
  iter <-  unlist(strsplit(temp, '_'))[5]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  DF_i <- read.csv(preds_iter[grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)])
  DF_i <- cbind(DF_i, preds_i[,c(setdiff(names(preds_i), names(DF_i)))])
  DF_i <- write.csv(DF_i, 
                    file = preds_iter[grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)],
                    row.names = FALSE)
}








# for RFs

preds_RFs <- list.files(path = 'out/', 
                        pattern = 'evalRF', 
                        full.names = TRUE)

for(i in 1:length(preds_RFs)){
  file_i <- preds_RFs[i]
  preds_RF_i <- read.csv(file_i)[7:14]
  temp <- tail(unlist(strsplit(file_i, '/')), n = 1)
  model <- unlist(strsplit(temp, '_'))[2]
  model <- paste0("RF_", model)
  iter <-  unlist(strsplit(temp, '_'))[3]
  iter <- as.numeric(unlist(strsplit(iter, '.csv')))
  iter_i <- grep(paste0('Eval_iter_', iter, '.csv'), x = preds_iter)
  preds_iter_i <- read.csv(preds_iter[iter_i])
  names(preds_RF_i) <- paste(model, x = names(preds_RF_i), sep = '_')
  preds_iter_i <- cbind(preds_iter_i, preds_RF_i)
  preds_iter_i <- write.csv(preds_iter_i, file = preds_iter[iter_i], row.names = FALSE)
}


