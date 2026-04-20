# Author: David E. Uribe-Rivera
# Does: create and export figures of the effects of contamination of PO data
# with different proportions of systematic surveys-derived records.
# Date: Dec 2025


library(ggplot2)
library(ggpubr)
library(terra)
library(RISDM)
library(tidyterra)
library(dplyr)

R.utils::sourceDirectory('R/')

##############################################################################
#   Figure 1: Maps of the study area and simulated species intensity
##############################################################################

# Figure will plot the simulation produced in iteration 1 (out of 100)
r <- rast('out/Sim_rast_1.tif')
plot(log(r$SppIntensity_1))

# read the controlPOdata
OpportunisticPO <- readRDS('out/control_cleanPO_1.rds')
OpportunisticPO_spat <- vect(as.data.frame(OpportunisticPO$PO), geom = c('x', 'y'), crs=r)
plot(OpportunisticPO_spat)

presencesPA <- as.data.frame(OpportunisticPO$PA)
absencesPA <- presencesPA[presencesPA$PA == 0,]
presencesPA <- presencesPA[presencesPA$PA == 1,]
absencesPA_spat <- vect(absencesPA, geom = c('x', 'y'), crs=r)
presencesPA_spat <- vect(presencesPA, geom = c('x', 'y'), crs=r)
points(absencesPA_spat, col = 'darkblue')
points(presencesPA_spat, col = 'red')
# reading the treatment of 50% contamination 
# (this means 2500 records come from the opportunistic PO data, 
# and the other 2500 records come from the systematic records ie PA or AA)
treatm1_50percent_1 <- readRDS('out/treatm1_50percent_1.rds')
PO_50perc_contam <- vect(as.data.frame(treatm1_50percent_1$PO), 
                         geom = c('x', 'y'),
                         crs = r)

# extract coordinates
xy_A <- crds(PO_50perc_contam, df = TRUE)
xy_B <- crds(OpportunisticPO_spat, df = TRUE)

# find rows in A not present in B
idx_not_in_B <- !paste(xy_A$x, xy_A$y) %in% paste(xy_B$x, xy_B$y)

# subset SpatVector
contamination <- PO_50perc_contam[idx_not_in_B, ]


# for plotting we will revert the UTM projection to map it in latlong units
r <- project(r, 'EPSG:4326')
contamination <- project(contamination, 'EPSG:4326')
OpportunisticPO_spat <- project(OpportunisticPO_spat, 'EPSG:4326')
PO_50perc_contam <- project(PO_50perc_contam, 'EPSG:4326')

# making the figure
pdf(file = 'out/Figure1_MapsContamination.pdf', width = 7, height = 7, onefile = TRUE)
par(mfrow = c(2, 2))
# plot panel A: Simulated species intensity
plot(log(r$SppIntensity_1))
north(crs = crs(r), len = 0.05, lab = "N", rotate = 0)

# plot panel B: Simulated sampling bias
plot(r$BiasIntensity_1)
north(crs = crs(r), len = 0.05, lab = "N", rotate = 0)

# plot panel C: Simulated (species intensity * sampling bias) and fully opportunistic PO records
plot(log(r$SppIntensity_1)+r$BiasIntensity_1)
north(crs = crs(r), len = 0.05, lab = "N", rotate = 0)
points(OpportunisticPO_spat, cex = 0.4, alpha = 0.7)

# plot panel D: Simulated species intensity and contaminated PO with 50% systematic records
plot(log(r$SppIntensity_1))
north(crs = crs(r), len = 0.05, lab = "N", rotate = 0)
points(contamination, col = 'red', cex = 0.4, alpha = 0.7)
points(OpportunisticPO_spat, col = 'darkblue', cex = 0.4, alpha = 0.7)

dev.off()


##############################################################################
#   Figure 2: coefficients estimates
##############################################################################

# boxplots of the posterior samples of the coefficient estimates for distribution
# linear and quadratic terms, and for linear bias-related covariates. 
# red dashed lines represent the parameter value used for simulating the data.

# First the complete set of coefficients for cSDMs without and with bias-related covariates, and for ISDMs
# these figures will be as supplementary material

Boxplots_Coefs_ISDMs <- boxplot_coeffs(coef_csv = 'out/coefs_master_ISDM.csv', # list of csvs with predictions and read as data.frame
                                       plot_file = 'out/FigureS2A_boxplots_coefs_RandomCont_ISDM.pdf')

Boxplots_Coefs_PPM_noBias <- boxplot_coeffs(coef_csv = 'out/coefs_master_PPM_noBias.csv', # list of csvs with predictions and read as data.frame
                                            plot_file = 'out/FigureS2B_boxplots_coefs_RandomCont_PPM_noBias.pdf')

Boxplots_Coefs_PPM_withBias <- boxplot_coeffs(coef_csv = 'out/coefs_master_PPM_withBias.csv', # list of csvs with predictions and read as data.frame
                                              plot_file = 'out/FigureS2C_boxplots_coefs_RandomCont_PPM_withBias.pdf')

# Now Figure 2, which will only display three of the eight distribution coefficients and both bias-coefficients
#put all six plots together into one multipanel plot

# bio1 (mean annual temp) - linear 
# set a common scale
Boxplots_Coefs_PPM_noBias$list_of_plots[[1]] <- Boxplots_Coefs_PPM_noBias$list_of_plots[[1]] + ylim(c(0.9, 2.8)) + 
  ggtitle("cSDM (PPM) - no bias-\nrelated covariates") + xlab("") + ylab("mean annual temp.\nlinear")
Boxplots_Coefs_PPM_withBias$list_of_plots[[1]] <- Boxplots_Coefs_PPM_withBias$list_of_plots[[1]] + ylim(c(0.9, 2.8)) +
  ggtitle("cSDM (PPM) - with bias-\nrelated covariates") + xlab("")
Boxplots_Coefs_ISDMs$list_of_plots[[1]] <- Boxplots_Coefs_ISDMs$list_of_plots[[1]] + ylim(c(0.9, 2.8)) +
  ggtitle("ISDM - explicitly \nmodelling SSB") + xlab("")

multi_bio1_linear <- ggarrange(Boxplots_Coefs_PPM_noBias$list_of_plots[[1]],
                        Boxplots_Coefs_PPM_withBias$list_of_plots[[1]],
                        Boxplots_Coefs_ISDMs$list_of_plots[[1]],
                        ncol = 3, #adjust plot space 
                        common.legend = T,
                        legend = 'none') #only legend in the last coef

# bio1 (mean annual temp) - quadratic 
# set a common scale
Boxplots_Coefs_PPM_noBias$list_of_plots[[2]] <- Boxplots_Coefs_PPM_noBias$list_of_plots[[2]] + ylim(c(-1.8, -0.75)) + 
  xlab("") + ylab("mean annual temp.\nquadratic") + ggtitle("")

Boxplots_Coefs_PPM_withBias$list_of_plots[[2]] <- Boxplots_Coefs_PPM_withBias$list_of_plots[[2]] + ylim(c(-1.8, -0.75)) + 
  xlab("") + ggtitle("")

Boxplots_Coefs_ISDMs$list_of_plots[[2]] <- Boxplots_Coefs_ISDMs$list_of_plots[[2]] + 
  ylim(c(-1.8, -0.75)) + xlab("") + ggtitle("")
  

multi_bio1_quadratic <- ggarrange(Boxplots_Coefs_PPM_noBias$list_of_plots[[2]],
                        Boxplots_Coefs_PPM_withBias$list_of_plots[[2]],
                        Boxplots_Coefs_ISDMs$list_of_plots[[2]],
                        ncol = 3, #adjust plot space 
                        common.legend = T,
                        legend = 'none') #only legend in the last coef


# bio12 (annual precipitation) - quadratic 
# set a common scale
Boxplots_Coefs_PPM_noBias$list_of_plots[[6]] <- Boxplots_Coefs_PPM_noBias$list_of_plots[[6]] + ylim(c(0, -1)) + 
  xlab("") + ylab("total precipitation \nquadratic") + ggtitle("")

Boxplots_Coefs_PPM_withBias$list_of_plots[[6]] <- Boxplots_Coefs_PPM_withBias$list_of_plots[[6]] + 
  ylim(c(0, -1)) + xlab("") + ggtitle("")

Boxplots_Coefs_ISDMs$list_of_plots[[6]] <- Boxplots_Coefs_ISDMs$list_of_plots[[6]] + 
  ylim(c(0, -1)) + xlab("") + ggtitle("")


multi_bio12_quadratic <- ggarrange(Boxplots_Coefs_PPM_noBias$list_of_plots[[6]],
                                  Boxplots_Coefs_PPM_withBias$list_of_plots[[6]],
                                  Boxplots_Coefs_ISDMs$list_of_plots[[6]],
                                  ncol = 3, #adjust plot space 
                                  common.legend = T,
                                  legend = 'none') #only legend in the last coef


# Bias-related covariates
# Human population density 
# set a common scale
Boxplots_Coefs_PPM_noBias$list_of_plots[[9]] <- Boxplots_Coefs_PPM_noBias$list_of_plots[[9]] + ylim(c(0.4, 1)) + 
  xlab("") + ylab("Human population \ndensity") + ggtitle("")

Boxplots_Coefs_PPM_withBias$list_of_plots[[9]] <- Boxplots_Coefs_PPM_withBias$list_of_plots[[9]] + 
  ylim(c(0.4, 1)) + xlab("") + ggtitle("")

Boxplots_Coefs_ISDMs$list_of_plots[[9]] <- Boxplots_Coefs_ISDMs$list_of_plots[[9]] + 
  ylim(c(0.4, 1)) + xlab("") + ggtitle("")

multi_HumanPop <- ggarrange(Boxplots_Coefs_PPM_noBias$list_of_plots[[9]],
                                   Boxplots_Coefs_PPM_withBias$list_of_plots[[9]],
                                   Boxplots_Coefs_ISDMs$list_of_plots[[9]],
                                   ncol = 3, #adjust plot space 
                                   common.legend = T,
                                   legend = 'none') #only legend in the last coef

# City accessibility 
Boxplots_Coefs_PPM_noBias$list_of_plots[[10]] <- Boxplots_Coefs_PPM_noBias$list_of_plots[[10]] + ylim(c(-1.5, -0.3)) + 
  xlab("") + ylab("City accessibility \n(sqrt)") + ggtitle("")

Boxplots_Coefs_PPM_withBias$list_of_plots[[10]] <- Boxplots_Coefs_PPM_withBias$list_of_plots[[10]] + 
  ylim(c(-1.5, -0.3)) + xlab("") + ggtitle("")

Boxplots_Coefs_ISDMs$list_of_plots[[10]] <- Boxplots_Coefs_ISDMs$list_of_plots[[10]] + 
  ylim(c(-1.5, -0.3)) + xlab("") + ggtitle("")

multi_CityAccess <- ggarrange(Boxplots_Coefs_PPM_noBias$list_of_plots[[10]],
                            Boxplots_Coefs_PPM_withBias$list_of_plots[[10]],
                            Boxplots_Coefs_ISDMs$list_of_plots[[10]],
                            ncol = 3, #adjust plot space 
                            common.legend = T,
                            legend = 'bottom') #only legend in the last coef


multi_plot <- ggarrange(multi_bio1_linear,
                        multi_bio1_quadratic,
                        multi_bio12_quadratic,
                        multi_HumanPop,
                        multi_CityAccess,
                        ncol = 1, common.legend = T,
                        legend = 'none') 

multi_plot<- annotate_figure(multi_plot,
                             left = grid::textGrob("Median coef. estimate", 
                                                   rot = 90, gp = grid::gpar(cex = 1.25)),
                             bottom = grid::textGrob("Percentage of survery-derived records in PO data", 
                                                     gp = grid::gpar(cex = 1.25)))

pdf(file = 'out/Figure2.pdf', 
    onefile=TRUE, width=9, height=11)
print(multi_plot)
dev.off()


# computing the rate of coef estimate failures (Credible intervals for coef 
# estimates do not cover the value used to simulate the data)

failure_rate_ISDMs_coefs <- list(
  Coefs = list(bio_1_lin = 1.7,  # bio_1_lin
               bio_1_quad = -1.4, # bio_1_quad
               bio_4_lin = 1.2, # bio_4_lin
               bio_4_quad = -1.5, # bio_4_quad
               bio_12_lin = -1.1, # bio_12_lin
               bio_12_quad = -0.5, # bio_12_quad
               bio_15_lin = -0.5, # bio_15_lin
               bio_15_quad = -1.7 ), #, # bio_15_quad
               # HumanPopulation = 0.9, # HumanPopulation
               # sqrtAccessibility = -1.2), # sqrtAccessibility
  coef_estimates = Boxplots_Coefs_ISDMs$coef_estimates,
  coef_fail_rate = array(data = 0, dim = c(8, 5, 2))
)

dimnames(failure_rate_ISDMs_coefs$coef_fail_rate)[[1]] <- names(failure_rate_ISDMs_coefs$Coefs)
dimnames(failure_rate_ISDMs_coefs$coef_fail_rate)[[2]] <- c('0%', 
                                                            '5%',
                                                            '10%', 
                                                            '25%',
                                                            '50%')
dimnames(failure_rate_ISDMs_coefs$coef_fail_rate)[[3]] <- c('no SRE', 
                                                            'with SRE')


 
for(i in names(failure_rate_ISDMs_coefs$Coefs)){
  basetemp <- dplyr::select(Boxplots_Coefs_ISDMs$coef_estimates, contains(i))
  j_0 <- 0
  for(j in c(0, 4:1)){
    j_0 <- j_0 + 1
    model <- paste0('m', j)
    for(k in c('', 'b')){
      freq <- freq1 <- freq2 <- 0
      if(k == ''){SRE <- 1}else{SRE <- 2} 
      model <- paste0(model, k)
      model2 <- paste0(model, '_')
      temp <- dplyr::select(basetemp, contains(model2))
      coef_id <- which(names(failure_rate_ISDMs_coefs$Coefs) == i)
      freq1 <- ifelse(length(table((temp[,1] > as.numeric(failure_rate_ISDMs_coefs$Coefs[coef_id])) )) == 2, 
                      table((temp[,1] > as.numeric(failure_rate_ISDMs_coefs$Coefs[coef_id])) )[2], 
                      table((temp[,1] > as.numeric(failure_rate_ISDMs_coefs$Coefs[coef_id])) )[1])
      freq2 <- ifelse(length(table((temp[,3] < as.numeric(failure_rate_ISDMs_coefs$Coefs[coef_id])) )) == 2, 
                      table((temp[,3] < as.numeric(failure_rate_ISDMs_coefs$Coefs[coef_id])) )[2], 
                      table((temp[,3] < as.numeric(failure_rate_ISDMs_coefs$Coefs[coef_id])) )[1])
      freq <- (freq1 + freq2)/2
      failure_rate_ISDMs_coefs$coef_fail_rate[i,j_0,SRE] <- as.numeric(freq)
    }
  }
}

write.csv(failure_rate_ISDMs_coefs$coef_fail_rate[,,1], file = "out/failure_rate_ISDMs_noSRE.csv")
write.csv(failure_rate_ISDMs_coefs$coef_fail_rate[,,2], file = "out/failure_rate_ISDMs_withSRE.csv")

# repeat for PPMs with no bias covariates

failure_rate_PPMs_noBias_coefs <- list(
  Coefs = list(bio_1_lin = 1.7,  # bio_1_lin
               bio_1_quad = -1.4, # bio_1_quad
               bio_4_lin = 1.2, # bio_4_lin
               bio_4_quad = -1.5, # bio_4_quad
               bio_12_lin = -1.1, # bio_12_lin
               bio_12_quad = -0.5, # bio_12_quad
               bio_15_lin = -0.5, # bio_15_lin
               bio_15_quad = -1.7 ), #, # bio_15_quad
  # HumanPopulation = 0.9, # HumanPopulation
  # sqrtAccessibility = -1.2), # sqrtAccessibility
  coef_estimates = Boxplots_Coefs_PPM_noBias$coef_estimates,
  coef_fail_rate = array(data = 0, dim = c(8, 5, 2))
)

dimnames(failure_rate_PPMs_noBias_coefs$coef_fail_rate)[[1]] <- names(failure_rate_PPMs_noBias_coefs$Coefs)
dimnames(failure_rate_PPMs_noBias_coefs$coef_fail_rate)[[2]] <- c('0%', 
                                                            '5%',
                                                            '10%', 
                                                            '25%',
                                                            '50%')
dimnames(failure_rate_PPMs_noBias_coefs$coef_fail_rate)[[3]] <- c('no SRE', 
                                                            'with SRE')



for(i in names(failure_rate_PPMs_noBias_coefs$Coefs)){
  basetemp <- dplyr::select(failure_rate_PPMs_noBias_coefs$coef_estimates, contains(i))
  j_0 <- 0
  for(j in c(0, 4:1)){
    j_0 <- j_0 + 1
    model <- paste0('m', j)
    for(k in c('', 'b')){
      freq <- freq1 <- freq2 <- 0
      if(k == ''){SRE <- 1}else{SRE <- 2} 
      model <- paste0(model, k)
      model2 <- paste0(model, '_')
      temp <- dplyr::select(basetemp, contains(model2))
      coef_id <- which(names(failure_rate_PPMs_noBias_coefs$Coefs) == i)
      freq1 <- ifelse(length(table((temp[,1] > as.numeric(failure_rate_PPMs_noBias_coefs$Coefs[coef_id])) )) == 2, 
                      table((temp[,1] > as.numeric(failure_rate_PPMs_noBias_coefs$Coefs[coef_id])) )[2], 
                      table((temp[,1] > as.numeric(failure_rate_PPMs_noBias_coefs$Coefs[coef_id])) )[1])
      freq2 <- ifelse(length(table((temp[,3] < as.numeric(failure_rate_PPMs_noBias_coefs$Coefs[coef_id])) )) == 2, 
                      table((temp[,3] < as.numeric(failure_rate_PPMs_noBias_coefs$Coefs[coef_id])) )[2], 
                      table((temp[,3] < as.numeric(failure_rate_PPMs_noBias_coefs$Coefs[coef_id])) )[1])
      freq <- (freq1 + freq2)/2
      failure_rate_PPMs_noBias_coefs$coef_fail_rate[i,j_0,SRE] <- as.numeric(freq)
    }
  }
}

write.csv(failure_rate_PPMs_noBias_coefs$coef_fail_rate[,,1], file = "out/failure_rate_PPMs_noBias_noSRE.csv")
write.csv(failure_rate_PPMs_noBias_coefs$coef_fail_rate[,,2], file = "out/failure_rate_PPMs_noBias_withSRE.csv")


# now for PPMs WITH bias covariates

failure_rate_PPMs_withBias_coefs <- list(
  Coefs = list(bio_1_lin = 1.7,  # bio_1_lin
               bio_1_quad = -1.4, # bio_1_quad
               bio_4_lin = 1.2, # bio_4_lin
               bio_4_quad = -1.5, # bio_4_quad
               bio_12_lin = -1.1, # bio_12_lin
               bio_12_quad = -0.5, # bio_12_quad
               bio_15_lin = -0.5, # bio_15_lin
               bio_15_quad = -1.7 ), #, # bio_15_quad
  # HumanPopulation = 0.9, # HumanPopulation
  # sqrtAccessibility = -1.2), # sqrtAccessibility
  coef_estimates = Boxplots_Coefs_PPM_withBias$coef_estimates,
  coef_fail_rate = array(data = 0, dim = c(8, 5, 2))
)

dimnames(failure_rate_PPMs_withBias_coefs$coef_fail_rate)[[1]] <- names(failure_rate_PPMs_withBias_coefs$Coefs)
dimnames(failure_rate_PPMs_withBias_coefs$coef_fail_rate)[[2]] <- c('0%', 
                                                                  '5%',
                                                                  '10%', 
                                                                  '25%',
                                                                  '50%')
dimnames(failure_rate_PPMs_withBias_coefs$coef_fail_rate)[[3]] <- c('no SRE', 
                                                                  'with SRE')



for(i in names(failure_rate_PPMs_withBias_coefs$Coefs)){
  basetemp <- dplyr::select(failure_rate_PPMs_withBias_coefs$coef_estimates, contains(i))
  j_0 <- 0
  for(j in c(0, 4:1)){
    j_0 <- j_0 + 1
    model <- paste0('m', j)
    for(k in c('', 'b')){
      freq <- freq1 <- freq2 <- 0
      if(k == ''){SRE <- 1}else{SRE <- 2} 
      model <- paste0(model, k)
      model2 <- paste0(model, '_')
      temp <- dplyr::select(basetemp, contains(model2))
      coef_id <- which(names(failure_rate_PPMs_withBias_coefs$Coefs) == i)
      freq1 <- ifelse(length(table((temp[,1] > as.numeric(failure_rate_PPMs_withBias_coefs$Coefs[coef_id])) )) == 2, 
                      table((temp[,1] > as.numeric(failure_rate_PPMs_withBias_coefs$Coefs[coef_id])) )[2], 
                      table((temp[,1] > as.numeric(failure_rate_PPMs_withBias_coefs$Coefs[coef_id])) )[1])
      freq2 <- ifelse(length(table((temp[,3] < as.numeric(failure_rate_PPMs_withBias_coefs$Coefs[coef_id])) )) == 2, 
                      table((temp[,3] < as.numeric(failure_rate_PPMs_withBias_coefs$Coefs[coef_id])) )[2], 
                      table((temp[,3] < as.numeric(failure_rate_PPMs_withBias_coefs$Coefs[coef_id])) )[1])
      freq <- (freq1 + freq2)/2
      failure_rate_PPMs_withBias_coefs$coef_fail_rate[i,j_0,SRE] <- as.numeric(freq)
    }
  }
}

write.csv(failure_rate_PPMs_withBias_coefs$coef_fail_rate[,,1], file = "out/failure_rate_PPMs_withBias_noSRE.csv")
write.csv(failure_rate_PPMs_withBias_coefs$coef_fail_rate[,,2], file = "out/failure_rate_PPMs_withBias_withSRE.csv")




##############################################################################
#   Figure 3: Model performance evaluations
##############################################################################

# boxplots for accuracy and discrimination ability of ISDMs & conventional SDMs


# i) accuracy measured by the absolute error of estimated occupancy probabilities
# (Error in median estimated occupancy probability: absolute errors = simulated - predicted values)

preds_CSVs <- list.files(path = 'out/', 
                         pattern = 'Eval_iter_', 
                         full.names = TRUE)

preds <- lapply(preds_CSVs, read.csv)

preds_ISDM <- lapply(preds, function(df) {
  select(df, contains("ISDM_"))
})

error_ISDM <- boxplot_absolute_errors(preds_list = preds_ISDM, # list of csvs with predictions and read as data.frame
                                        plot_identifier = 'OccProb',
                                        model_type = 'ISDM',
                                        plot_file = "out/Fig3_RandomCont_ISDM_OccProb_Error.pdf")
error_ISDM <- error_ISDM + ylim(c(0, 0.13)) + 
  ylab("") + ggtitle("ISDM - explicitly \nmodelling SSB") + xlab("") +
  theme(plot.title = element_text(size = 10))


preds_PPM_noBias <- lapply(preds, function(df) {
  select(df, contains("PPM_noBias"))
})
error_PPM_noBias <- boxplot_absolute_errors(preds_list = preds_PPM_noBias, # list of csvs with predictions and read as data.frame
                        plot_identifier = 'OccProb',
                        model_type = 'ISDM',
                        plot_file = "out/Fig3_RandomCont_PPM_noBias_OccProb_Error.pdf")
error_PPM_noBias <- error_PPM_noBias + ylim(c(0, 0.13)) +
  ggtitle("cSDM (PPM) - no bias-\nrelated covariates") + xlab("") +
  theme(plot.title = element_text(size = 10))


preds_PPM_withBias <- lapply(preds, function(df) {
  select(df, contains("PPM_withBias"))
})
error_PPM_withBias <- boxplot_absolute_errors(preds_list = preds_PPM_withBias, # list of csvs with predictions and read as data.frame
                                            plot_identifier = 'OccProb',
                                            model_type = 'ISDM',
                                            plot_file = "out/Fig3_RandomCont_PPM_withBias_OccProb_Error.pdf")
error_PPM_withBias <- error_PPM_withBias + ylim(c(0, 0.13)) + 
 ggtitle("cSDM (PPM) - with bias-\nrelated covariates") + ylab("") + xlab("") +
  theme(plot.title = element_text(size = 10))






multi_plot<- ggarrange(error_PPM_noBias,error_PPM_withBias, error_ISDM,
                       ncol = 3, common.legend = T,
                       legend = 'bottom') 


# preds_RF <- lapply(preds, function(df) {
# select(df, contains("RF"))
# })
# 
# # now adding AUC-ROC evaluations for both models
# AUCROC_plots <- boxplot_AUCs(results_csv = 'out/results_master.csv', 
#              pref_sampling = FALSE,
#              plot_file = 'out/boxplots_AUCs_RandomCont.pdf')
# 
# 
# multi_plot<- ggarrange(AUCROC_plots,multi_plot,
#                        ncol = 1, common.legend = T,
#                        legend = 'top')
# 
multi_plot<- ggpubr::annotate_figure(multi_plot,
                                     # left = grid::textGrob("Median absolute error",
                                                           # rot = 90, gp = grid::gpar(cex = 1.25)),
                                     bottom = grid::textGrob("Percentage of survery-derived records in PO data",
                                                             gp = grid::gpar(cex = 1)))

pdf(file = 'out/Figure3.pdf', width = 8, height = 5, onefile = TRUE)
multi_plot
dev.off()

