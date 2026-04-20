# testing RISDM behaviour on simulated data when PO data is mixed with survey-derived records

# for HPC
.libPaths("/home/uri003/EffectsofsystenaticRecsinPOdata/lib/")

#Read the command line arguments
command_args <- commandArgs(trailingOnly = TRUE)

# Define command line arguments as a variable
i <- command_args[1]
i <- as.numeric(i)
print(" # set.seed = ")
print(i)

# 1. 30 replicates over this script, where each replicate is a different
#   realisation, the different between rounds is driven by a spatial random effect.

library(RISDM)
library(terra)
library(predicts)
library(enmSdmX)
library(randomForest)

R.utils::sourceDirectory('R/')

covars <- rast('data/Tas_env.tif') # at ~0.88 Km resolution
names(covars)[5:6] <- c('HumanPopulation', 'sqrtAccessibility')
# making orthogonal polynomies
covars <- c(make_orthogonal_poly(covars[[1:4]]), covars[[5:6]])
names(covars)
covars <- scale(covars)
plot(covars)

distForm <- formula(~ -1 + bio_1_lin + bio_1_quad + 
                      bio_4_lin + bio_4_quad +
                      bio_12_lin + bio_12_quad +
                      bio_15_lin + bio_15_quad)

biasForm <- formula(~ 1 + HumanPopulation + sqrtAccessibility)

# 2. use risdm to simulate individuals across space based on realisation of
#   point process
dat_withRE <- simulateData.isdm(
  pop.size = 800000,
  n.PO = 5000,
  n.PA = 16000,
  n.AA = 16000,
  n.DC = 100,
  distForm = distForm,
  biasForm = biasForm,
  distCoefs = c(1.7,  # bio_1_lin
                -1.4, # bio_1_quad
                1.2, # bio_4_lin
                -1.5, # bio_4_quad
                -1.1, # bio_12_lin
                -0.5, # bio_12_quad
                -0.5, # bio_15_lin
                -1.7), # bio_15_quad
  biasCoefs = c(-17, # intercept
                0.9, # HumanPopulation
                -1.2), # sqrtAccessibility
  
  DC.pi = matrix(c(0.8,0.76, 0.9,0.85, 0.82,0.87),
                 nrow = 3, ncol = 2, byrow = TRUE),
  transect.size = 0.25, #a proportion of cell size.
  covarBrick = covars[[c('bio_1_lin', 'bio_1_quad', 'bio_4_lin', 'bio_4_quad', 
                         'bio_12_lin',  'bio_12_quad', 'bio_15_lin', 'bio_15_quad',
                         'HumanPopulation', 'sqrtAccessibility')]],
  control = list(
    # addRandom = FALSE,
    range = 25 * max(terra::res(covars)),  # default in RISDM original simulation function
    set.random.seed = TRUE,
    random.seed = i,
    sd = 0.5,
    # useMBHdesign = TRUE,
    # exact.n.PO = TRUE,
    doPlot = FALSE)
)

plot(dat_withRE$covarBrick$LinPred)
plot(dat_withRE$covarBrick$Intensity)
# plot(dat_withRE$covarBrick$REff)
# plot(dat_withRE$covarBrick$biasIntensity)
#check visually the simulated relative abundance of koalas
round(sum(dat_withRE$covarBrick$Intensity[], na.rm = TRUE))
# 800,000

# saving the SpatRaster for each iteration with only the relevant layers:
# virtual species intensity, sampling bias intensity and the spatial random effect

r <- c(dat_withRE$covarBrick$Intensity, 
       (dat_withRE$covarBrick$LinPred-dat_withRE$covarBrick$biasLinPred),
       dat_withRE$covarBrick$REff)

names(r) <- c('SppIntensity', 'BiasIntensity', 'REff')
names(r) <- paste0(names(r), '_', i)
writeRaster(r, file = paste0('out/Sim_rast_', i, '.tif'),
            overwrite = TRUE)


# 3. uses simulated data to train 9 models with different degrees of survey
#   derived records into presence only. 4 levels of contamination then control.
#   2 types of contamination (random v blocked)
dat_random_structured <- dat_withRE[2:3]
dat_random_control_cleanPO <- list(
  PO = dat_withRE$PO,
  PA = head(dat_random_structured$PA, n = 500),
  AA = head(dat_random_structured$AA, n = 500)
)
str(dat_random_control_cleanPO)
saveRDS(dat_random_control_cleanPO, file = paste0('out/control_cleanPO_', i, '.rds'))


dat_random_structured_notused <- list(
  PA = dat_random_structured$PA[501:15500,],
  AA = tail(dat_random_structured$AA, n = 15500)
)
str(dat_random_structured_notused)

# out-of sample evaluation dataset
PA_eval <- tail(dat_withRE$PA, n = 500)
saveRDS(PA_eval, file = paste0("out/PA_eval_", i, ".rds"))

# treatment 1: 50% PO data coming from surveys
dat_random_treatm1_50percent <- create_treatment_dataset(
  percent_contaminated = 50, blocked = FALSE, iter = i)

# treatment 2: 25% PO data coming from surveys
dat_random_treatm2_25percent <- create_treatment_dataset(
  percent_contaminated = 25, blocked = FALSE, iter = i)

# treatment 3: 10% PO data coming from surveys
dat_random_treatm3_10percent <- create_treatment_dataset(
  percent_contaminated = 10, blocked = FALSE, iter = i)

# treatment 4: 5% PO data coming from surveys
dat_random_treatm4_5percent <- create_treatment_dataset(
  percent_contaminated = 5, blocked = FALSE, iter = i)





############################################################################
# Running the 10 models and saving the parameter estimates for each one


# mesh parameters
mesh.max.n <- c(8000, 1000)

meshy <- makeMesh(ras = dat_withRE$covarBrick$REff,
                  max.n = mesh.max.n,
                  dep.range = 10000,
                  expans.mult = 5,
                  offset = 100000,
                  doPlot = FALSE)



# m0: Control model (no contamination with survey data)
loop_spatial_re(
  dataset = dat_random_control_cleanPO,
  model_name = 'm0',
  Env_raster = covars,
  distForm = distForm,
  iter = i,
  mesh = meshy,
  pa_eval = PA_eval
)

# m1: 50% survey data in PO model (no contamination with survey data)
loop_spatial_re(
  dataset = dat_random_treatm1_50percent,
  model_name = 'm1',
  distForm = distForm,
  Env_raster = covars,
  iter = i,
  mesh = meshy,
  pa_eval = PA_eval
)

# m2 = 25% survey data in the PO dataset
loop_spatial_re(
  dataset = dat_random_treatm2_25percent,
  model_name = 'm2',
  distForm = distForm,
  Env_raster = covars,
  iter = i,
  mesh = meshy,
  pa_eval = PA_eval
)

# m3 = 10% survey data in the PO dataset
loop_spatial_re(
  dataset = dat_random_treatm3_10percent,
  model_name = 'm3',
  distForm = distForm,
  Env_raster = covars,
  iter = i,
  mesh = meshy,
  pa_eval = PA_eval
)

# m4: 5% survey data in PO model (no contamination with survey data)
loop_spatial_re(
  dataset = dat_random_treatm4_5percent,
  model_name = 'm4',
  distForm = distForm,
  Env_raster = covars,
  iter = i,
  mesh = meshy,
  pa_eval = PA_eval
)