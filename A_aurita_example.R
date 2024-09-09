# A script illustrating the use of the JellySim, DMat, and SimPlot functions in predicting spatial
# and temporal patterns in the movements of Aurelia aurita individuals within the Baltic Sea

# Primary Author: James Cant
# Date Last modified: Dec 2023
# ------------------------------------------------------------------------------------

##############################################
# STEP 1: Load required files
##############################################

# Define file path
filePath <- 'file_directory_path/'

# Load required functions
source(paste0(filePath, 'JellySim.R'))
source(paste0(filePath, 'DMat.R'))
source(paste0(filePath, 'SimPlot.R'))
source(paste0(filePath, 'ParSens.R'))

# Load demographic parameters
demogpars <- read.csv(paste0(filePath, 'Demographic parameters.csv'), row.names = 1)
demogpars <- with(demogpars, setNames(x, rownames(demogpars))) # ensure appropriate row naming

# Load drifting details
load(paste0(filePath, 'ParticleTracks2021.RData'))


##############################################
# STEP 2: Define model parameters
##############################################

# define user specific model parameters
# temporal resolutions of drifting simulations
nd = 30
# simulation spatial extents
xx = 31
xn = 3
yx = 66
yn = 53
# selected release locations
sites = c(598,240,954,597) # corresponds with the selection of sites in the central Baltic, the eastern Baltic, the Gulf of Finland and the Skagerrak & Kattegat region. 
# release timings
months = c(2:5) # corresponds with February through to May inclusive
# Demographic matrix dimensions
m = 200
# Resampling iterations
z = 10
# Multicore processing (True or False)
pl = TRUE


##############################################
# STEP 3: Run simulation
##############################################

JellyData <- JellySim(# demographic details
                      pars = demogpars, 
                      m = m,
                      # drifting details
                      driftData = SINMOD_data, 
                      n_days = nd,
                      # spatial details
                      xmx = xx, xmn = xn, ymx = yx, ymn = yn,
                      rel_location = sites,
                      # temporal details
                      rel_months = months,
                      # Resampling details
                      zmax = z,
                      # request for multicore processing
                      parallel = pl)


##############################################
# STEP 4: Visualize Outputs
##############################################

JellyPlots <- SimPlot(# main forecasts
                       meanRast = JellyData$mean, 
                       # plots showing forecast confidence
                       confRast = JellyData$conf,
                       # Highlight initial release sites
                       sites = JellyData$site,
                       # Spatial extent of plots
                       xmx = xx, xmn = xn, ymx = yx, ymn = yn)

# Mean forecasts
JellyPlots$mean

# Forecast confidence
JellyPlots$conf


##############################################
# STEP 5: Quantify parameter sensitivities
##############################################

# Resupply the same details to the ParSens arguments as provided to JellySim above, 
# less a value for zmax and selecting only one of the four release locations
JellySens <- ParSens(pars = demogpars, m = m, driftData = SINMOD_data, n_days = nd,
                     xmx = xx, xmn = xn, ymx = yx, ymn = yn, rel_location = sites, rel_months = months, parallel = pl,
                     params = c('P.dens.temp', 'M.surv.temp', 'Growth.slope', 'E.surv.temp', 'No.ephyra.temp'))
                            # Polyp Density, Medusae survival, Medusae growth, Ephyra survival, Number of Ephyra produced.
# Sensitivities
JellySens$plot

# ------------------------------------------------ End of Code ----------------------------------