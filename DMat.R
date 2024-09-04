# This script outlines the DMat function used within the GoJelly predictive tool for forecasting the spatial dynamics of jellyfish blooms.
# The primary purpose of this function is to generate periodic population matrix models using a series of pre-estimated demographic variables and environmental covariates.

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
# -----------------------------------------------------------------------------

#----------------------------------------------------
# STEP 1: DEFINE PRIMARY FUNCTION
#----------------------------------------------------

# DMat serves as a demographic black box. It is provided with a series of abiotic readings and demographic parameters describing the various survival, growth and reproductive patterns
# underlying the dynamics of the population of interest, and how these patterns are influenced by abiotic conditions. This is a flexible component of the GoJelly predictive tool and can be modified 
# depending on the unique dynamics of the focal population. Crucially, all that needs to be maintained is that function must return a list comprising a sequence of survival and growth matrices (matU),
# clonality matrices (matC), reproductive matrices (matF), combined population matrices (matB, equal to matU+matC+matF), a measure of initial population size, 
# the widths of any size bins used (continuous matrix models), and an indexing factor needed to align monthly demographic rates and abiotic covariates.

# Function for generating combined discrete-continuous periodic population matrix
DMat <- function(m, n_month, pars, Temp, Sal, rel_months, ephyra_omit) {
   
  # Load package dependencies
  packages <- c("purrr", "gamlss", "gamlss.dist")
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    invisible(install.packages(packages[!installed_packages]))
  }
  invisible(lapply(packages, library, character.only = TRUE))
  
  ### This function needs vectors describing the temperatures & salinities experienced: 1. by polyps at point of ephyra release, and 2. the subsequent drifting medusae.
  # Define a value to keep indexing consisted across simulation runs for when the number of periods being simulated is less than the total temporal expanse of the simulation
  index <- length(Temp[is.na(Temp)])+1
  # shift temperature vector to ensure it matches formatting
  Temp <- data.table::shift(Temp, n = index-1, fill = NA)
  Sal <- data.table::shift(Sal, n = index-1, fill = NA)
  # And remove error entries
  Temp[Temp == -100] <- NA
  Sal[Sal == -100] <- NA
  # Define a vector for vectorising loop functions
  index_vec <- 1:(n_month) # each month in the simulation corresponds with 30 days.
  
  ### 1. Start by calculating the initial polyp density (per km2) for this simulation
  densPolyp <- P_dens(pars, Temp = Temp[index], Sal = Sal[index])
  
  ### 2. Next define the vital rates/variables that are not temperature nor time dependent
  # larval stasis
  Plan_stasis <- as.numeric(rep(pars["Plan_stasis"], length = n_month))
  # minimum medusae size
  L <- pars["min.size"]*0.9 # smallest observed size minus 10% to avoid eviction.
  
  ### 3. define periodic estimates for discrete transition values/population variables
  Plan_Polyp <- replicate(n_month, Plan_settle(pars))
  Polyp_Polyp <- sapply(index_vec, P_stasis, pars = pars, Temp = Temp, Sal = Sal)
  Ephyra_Ephyra <- sapply(index_vec, E_stasis, pars = pars, Temp = Temp)
  Polyp_Ephyra <- sapply(index_vec, P_Ephyra, pars = pars, Psurv = Polyp_Polyp, Temp = Temp)
  Polyp_bud <- sapply(index_vec, P_bud, pars = pars, Psurv = Polyp_Polyp, Temp = Temp)
  Ripe_size <- sapply(index_vec, Ripe_z, pars = pars, Temp = Temp)
  # define medusae survival - although this is associated with a continuous stage the vital rate itself is discrete
  Med_Med <- sapply(index_vec, s_z, pars = pars, Temp = Temp, Sal = Sal)
  # similarly max medusae size is also a discrete periodic entity associated with the continuous stage
  max.size.mean <- map_dbl(Temp, ~pars["max.size.int"] + (pars["max.size.temp"] * .x))
  U.store <- suppressWarnings(map_dbl(max.size.mean, ~exp(rnorm(1, .x, pars["max.size.sd"])) * 1.1)) # largest possible size plus 10% and correct for log scale used

  rm(max.size.mean)
  # remove NaNs
  U.store[is.nan(U.store)] <- NA
  
  ### 4. Parameterise the periodic mesh points which will be used by the continuous stage of the model to determine size for each periodic matrix
  # number of bins
  h <- map_dbl(U.store, ~(.x - L)/m)
  # meshpoints 
  meshpts <- map(h, ~L + ((1:m) - 1/2) * .x)
  
  ### 5. Parameterise periodic estimates for continuous transition values/population variables
  Plan_Polyp_Medusa_store <- lapply(index_vec, Ephyra_Medusa, Sizet1 = meshpts, pars = pars, Esurv = rep(0, length.out = length(index_vec)), Temp = Temp)
  Ephyra_Medusa_store <- lapply(index_vec, Ephyra_Medusa, Sizet1 = meshpts, pars = pars, Esurv = Ephyra_Ephyra, Temp = Temp)
  # create storage for remaining periodic values
  Pz1z_store <- Fz1z_store <- list()
  # loop through each period, extracting the associated abiotic conditions and using them to estimate periodic continuous transition matrices/vectors
  for(ii in index:(n_month)){
    Pz1z_store[[ii]] <- outer(meshpts[[ii]], meshpts[[ii]], Pz1z, pars = pars, survM = Med_Med[ii])
    Fz1z_store[[ii]] <- as.numeric(sapply(meshpts[[ii]], Fz1z, pars = pars, Ripe_size = Ripe_size[ii], Temp = Temp[ii], Sal = Sal[ii], survM = Med_Med[ii]))
  }
  
  ### 6. For each monthly model it is necessary to inhibit some vital rates (even if abiotic conditions allow them) to reflect temporal delays across the annual life-cycle of Aurelia
  # i.e. not all processes happen every month!
  # data in the biobank and Goldstein and Steiner 2020 indicates that larval settlement and persistence within the water column 
  # is seasonal and so needs to be prevented in certain months.
  # Based on past observations Larval settlement only occurs between July and December (Periods/Months 7-12)
  Plan_Polyp[1:6] <- 0
  # Similarly larvae only persist within the water column between August and Dec (periods 8-12)
  Plan_stasis[1:7] <- 0
  # the remaining functions either are variable year round (following abiotic regulation) or appear seasonal but are strongly inhibited by thermal regimes
  # so I do not need to artificially constrain them here.
  
  ### 7. Define each monthly population matrix (Bn)
  # Each period matrix will comprise combined discrete to continuous population models each consisting of four separate elements that will be stitched together
  # 1. A discrete 4*4 matrix detailing the transitions of planulae larvae, polyps, and ephyra.
  #    IMPORTANT: This matrix contains two Polyp stages. One is the seed bank polyps which continually supply ephyra, 
  #    the second is the polyps produced by drifting medusae to prevent them reseeding back into initial polyp population as in reality the larvae wouldn't travel back
  # 2. A discrete to continuous transition detailing the development of ephyra to medusae
  # 3. A continuous element describing the dynamics of medusae
  # 4. A continuous to discrete transition detailing the production of planulae larvae by medusae
  # First build the various discrete and continuous components of each periodic model before stitching them together.
  
  # Condense together all vital rate transitions associated with each month.
  ### Survival
  # Stitch together survival & growth transitions
  # Panel 1 (Discrete matrix element)
  U1 <- lapply(index_vec, mk_mat1, mat = "U", dims = 4, index = index, Plan_stasis, Plan_Polyp, Polyp_Polyp, Ephyra_Ephyra, Polyp_bud, Polyp_Ephyra)
  # Panel 2 (Discrete to continuous)
  U2 <- lapply(index_vec, mk_mat2, m=m, mat = "U", dims = 4, index = index, Plan_Polyp_Medusa_store, Ephyra_Medusa_store)
  # Panel 3 (Continuous Medusae survival and growth)                     
  U3 <- lapply(index_vec, mk_mat3, mat = "U", index = index, h, m = m, Pz1z_store)
  # Panel 4 (Continuous to Discrete)
  U4 <- lapply(index_vec, mk_mat4, m=m, mat = "U", index = index, dims = 4, h, Fz1z_store)
  
  # and repeat for clonality and reproduction
  # matC - Clonality
  C1 <- lapply(index_vec, mk_mat1, mat = "C", dims = 4, index = index, Plan_stasis, Plan_Polyp, Polyp_Polyp, Ephyra_Ephyra, Polyp_bud, Polyp_Ephyra)
  C2 <- lapply(index_vec, mk_mat2, m=m, mat = "C", dims = 4, index = index, Plan_Polyp_Medusa_store, Ephyra_Medusa_store)
  C3 <- lapply(index_vec, mk_mat3, mat = "C", index = index, h, m = m, Pz1z_store)
  C4 <- lapply(index_vec, mk_mat4, m=m, mat = "C", dims = 4, index = index, h, Fz1z_store)
  
  # matF - Fecundity
  F1 <- lapply(index_vec, mk_mat1, mat = "F", dims = 4, index = index, Plan_stasis, Plan_Polyp, Polyp_Polyp, Ephyra_Ephyra, Polyp_bud, Polyp_Ephyra)
  F2 <- lapply(index_vec, mk_mat2, m=m, mat = "F", dims = 4, index = index, Plan_Polyp_Medusa_store, Ephyra_Medusa_store)
  F3 <- lapply(index_vec, mk_mat3, mat = "F", index = index, h, m = m, Pz1z_store)
  F4 <- lapply(index_vec, mk_mat4, m=m, mat = "F", dims = 4, index = index, h, Fz1z_store)
  
  # Now stitch together the various sub matrix panels and condense into the monthly matrices
  # Survival
  matU <- sapply(index_vec, build_array,  mat = "U", U1 = U1, U2 = U2, U3 = U3, U4 = U4, simplify = 'array')
  # Clonality
  matC <- sapply(index_vec, build_array,  mat = "C", C1 = C1, C2 = C2, C3 = C3, C4 = C4, simplify = 'array')
  # Fecundity
  matF <- sapply(index_vec, build_array,  mat = "F", F1 = F1, F2 = F2, F3 = F3, F4 = F4, simplify = 'array')
  # And finally create the monthly matrix (Bring together survival, clonality, and growth)
  Bmat <- sapply(index_vec, mk_matB, matU, matC, matF, simplify = 'array')
  
  # However, before these matrices can be return they needed to be modified to adhere to the temporal constraints associated with ephyra release timings (i.e. not all transitions are permisable at all times)
  # Firstly, all ephyra production transitions (element 4,2) in all matrices corresponding with periods outside the user defined release window need to be set to zero
  Bmat[4,2,ephyra_omit] <- 0
  # Following initial Ephyra release the simulation then splits into two components:
  # The first component is only interested in the transitions and survival of produced Ephyra/Medusae.
  # the second is interested in how the dynamics of the polyps affects them for subsequent months.
  # Following initial ephyra release, all subsequent matrices need to have all transitions prior to ephyra production removed.
  # This can be achieved by removing matrix elements from column 1:3
  # First the function needs to define some indexing values for subsetting ephrya transitions (primarily to prevent subscript issues)
  ephyra.index1 <- min(rel_months)+1
  if(ephyra.index1 > n_month){ephyra.index1 = n_month}
  ephyra.index2 <- min(rel_months)+3
  if(ephyra.index2 > n_month){ephyra.index2 = n_month}
  ephyra.index3 <- index+1
  if(ephyra.index3 > n_month){ephyra.index3 = n_month}
  ephyra.index4 <- index+3
  if(ephyra.index4 > n_month){ephyra.index4 = n_month}
  ephyra.index5 <- max(rel_months) + 3
  if(ephyra.index5 > n_month){ephyra.index5 = n_month}
  # now apply inhibition.
  if(index < min(rel_months)){Bmat[,1:3, ephyra.index1:n_month] <- 0}
  if(index %in% rel_months){Bmat[,1:3, ephyra.index3:n_month] <- 0}
  if(index > max(rel_months)){Bmat[,1:3, ephyra.index1:n_month] <- 0}
  # Additionally, following their initial release Ephyra are only believed to persist in the water column for 2-3 months
  # Therefore this constraint needs to be implemented within the relevant matrices by artificially setting ephyra survival (element 4,4), and medusae production (elements 5:204,4) to zero in matrices corresponding with
  # periods 2-3 months after initial release.
  if(index < min(rel_months)){Bmat[4:(m+4),4, ephyra.index2:n_month] <- 0}
  if(index %in% rel_months){Bmat[4:(m+4),4, ephyra.index4:n_month] <- 0}
  if(index > max(rel_months)){Bmat[4:(m+4),4, ephyra.index5:n_month] <- 0}
  # Also medusae production is not possible in the first month of ephyra production and so needs inhibiting
  if(index < min(rel_months)){Bmat[5:(m+4),4, min(rel_months)] <- 0}
  if(index %in% rel_months){Bmat[5:(m+4),4, index] <- 0}
  if(index > max(rel_months)){Bmat[5:(m+4),4,  min(rel_months)] <- 0}
  
  # and return outputs
  return(list(matB = Bmat, matU = matU, matC = matC, matF = matF, DensityStart = densPolyp, h = h, Adjust = index))
}

#----------------------------------------------------
# STEP 2: DEFINE INTERNAL FUNCTIONS 
#----------------------------------------------------

# These functions are called internally within the DMat function above and correspond with each of the vital rate patterns that require paramterisation in order to 
# generate the required matrix population model. These functions can be modified and/or removed depending on the life-cycle of the focal population.

# 1. Define the functions for creating matrix models
# generate transition probabilities
define_probs <- function(x_bar, x_sd) {
  # randomly generated value from modeled distribution
  linearp <- rnorm(1, x_bar, x_sd)
  # reverse logit link
  p <- 1/(1+exp(-linearp))
  # quick fix to prevent impossible probabilities
  if(p > 1) {p = 1}
  if(p < 0) {p = 0}
  # extract probability
  return(p)
}

# convert gamlss Zero adjusted gamma model outputs into individual counts. -----------------------------------------------------------
# Zero-inflated gamma models are designed around the changing probability of zero entries.
convert_gamlss <- function(mu, sigma, nu) {
  # convert nu parameter into probability
  nu_p <- 1/(1+exp(-nu))
  # randomly generated value from modeled distribution
  count_n <- rZAGA(1, mu = exp(mu), sigma = exp(sigma), nu = nu_p)
  # extract count
  return(count_n)
}

# 2. Define function for estimating initial Polyp density -----------------------------------------------------------
P_dens <- function(pars, Temp, Sal){
  if(is.na(Temp)){
    densP <- NA
  } else {
    # polyp density as a function of temperature
    mean.dens <- pars["P.dens.int"] + (Temp * pars["P.dens.temp"]) +
      (Sal * pars["P.dens.sal"]) + (Temp * Sal * pars["P.dens.tempsal"])
    lineardensP <- rnorm(1, mean.dens, pars["P.dens.sd"])
    # and adjust for gamma transformation to give polyp density per cm2
    densP <- exp(lineardensP)
  }
  # and return polyp density
  return(densP)
}

# 3. Define vital rate functions associated with discrete stage transitions -----------------------------------------------------------
# Planulae settlement probability
Plan_settle <- function(pars){
  settle <- define_probs(pars["Settle.mean"], pars["Settle.sd"])
  return(settle)
}

# Polyp survival probability -----------------------------------------------------------
P_stasis <- function(ii, pars, Temp, Sal){
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    stasis <- NA
  } else {
    # survival as a function of temperature
    surv.mean <- pars["P.surv.int"] + (Temp[ii] * pars["P.surv.temp"])
    surv <- define_probs(surv.mean, pars["P.surv.sd"])
    # mortality adjustment in response to salinity
    mort.mean <- pars["P.mort.int"] + (Sal[ii] * pars["P.mort.sal1"]) + ((Sal[ii]^2) * pars["P.mort.sal2"])
    add_mort <- define_probs(mort.mean, pars["P.mort.sd"])
    # combine together
    stasis <- surv - add_mort
    # quick fix to prevent impossible survival probabilities
    if(stasis > 1) {stasis = 1}
    if(stasis < 0) {stasis = 0}
  }
  # return desired output
  return(stasis)
}

# Ephyra survival probability -----------------------------------------------------------
E_stasis <- function(ii, pars, Temp){
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    surv <- NA
  } else {
    # survival as a function of temperature
    surv.mean <- pars["E.surv.int"] + (Temp[ii] * pars["E.surv.temp"])
    surv <- define_probs(surv.mean, pars["E.surv.sd"])
    # quick fix to prevent impossible survival probabilities
    if(surv > 1) {surv = 1}
    if(surv < 0) {surv = 0}
  }
  return(surv)
}

# Ephyra production from Polyps (Strobilation) -----------------------------------------------------------
P_Ephyra <- function(ii, pars, Psurv, Temp){ # this will be the survival already generated for Polyps (so they match across transitions)
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    P_ephyra <- NA
  } else {
    # Probability of strobilation as a function of temperature
    strob.mean <- pars["P.strob.int"] + (Temp[ii] * pars["P.strob.temp"])
    P.strob <- define_probs(strob.mean, pars["P.strob.sd"])
    # Number of Ephyra produced as a function of temperature
    No.ephyra.mean <- pars["No.ephyra.int"] + (Temp[ii] * pars["No.ephyra.temp"])
    No.ephyra <- rnorm(1, No.ephyra.mean, pars["No.ephyra.sd"])

    # Combine ephyra production with the probability of strobilation and Polyp survival
    P_ephyra <- Psurv[ii] * P.strob * No.ephyra
    # just a quick fix to prevent negative production
    if(P_ephyra < 0) {P_ephyra = 0}
  }
  # and return
  return( P_ephyra )
}

# Polyp asexual reproduction -----------------------------------------------------------
P_bud <- function(ii, pars, Psurv, Temp){ # this will be the survival already generated for Polyps (again so they match across transitions)
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    bud <- NA 
    return( bud )
  } else {
    # Because this function is based on zero-inflated regression models it is not possible to remove variability here as the values output from the model
    # are dependent on the probability of zero entries.
    # Number of lateral buds produced as a function of temperature
    lat.mu.mean <- pars["No.lat.P.mu"] + (Temp[ii] * pars["No.lat.P.mu.temp"])
    No.lat <- convert_gamlss(lat.mu.mean, pars["No.lat.P.sigma"], pars["No.lat.P.nu"])
    # Number of stolons produced as a function of temperature
    stol.mu.mean <- pars["No.stol.P.mu"] + (Temp[ii] * pars["No.stol.P.mu.temp"])
    No.stol <- convert_gamlss(stol.mu.mean, pars["No.stol.P.sigma"], pars["No.stol.P.nu"])
    # Number of stolon buds produced as a function of temperature
    both.mu.mean <- pars["No.both.P.mu"] + (Temp[ii] * pars["No.both.P.mu.temp"])
    No.both <- convert_gamlss(both.mu.mean, pars["No.both.P.sigma"], pars["No.both.P.nu"])
    # combine together with polyp survival and likelihood of asexual reproduction
    return( as.numeric(Psurv[ii] * pars["P.Asexual"] * (No.lat + No.stol + No.both)) )
  }
}

# Size at which medusae achieve sexual maturity -----------------------------------------------------------
Ripe_z <- function(ii, pars, Temp) {
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    Ripe_size <- NA
  } else {
    mean <- pars["Ripe.size.int"] + (pars["Ripe.size.temp"] * Temp[ii]) 
    sd <- pars["Ripe.size.sd"]
    # randomly generated value from modeled distribution
    Ripe_size <- exp(rnorm(1, mean, sd)) # correct for the log scale used.
  }
  # return output
  return(Ripe_size)
}

# Medusae survival probability (there is no data for this to be related to medusae size but it is being -----------------------------------------------------------
# constrained by known physiological thresholds of A. aurita)
s_z <- function(ii, pars, Temp, Sal){
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    p <- NA
  } else {
    # survival as a function of temperature
    mean.surv <- pars["M.surv.int"] + (pars["M.surv.temp"] * Temp[ii])
    p <- define_probs(mean.surv, pars["M.surv.sd"])
    # quick fix to prevent impossible probabilities
    if(p > 1) {p = 1}
    if(p < 0) {p = 0}
    # and a little fix to ensure survival doesn't occur below certain temperatures and salinities
    if (Temp[ii] < pars["surv.min.temp"] || Sal[ii] < pars["surv.min.sal"])  { p = 0 }
  }
  # return output
  return(p) 
}

# 4. Define functions associated with continuous state transitions. -----------------------------------------------------------
# Medusae production from Ephyra (metamorphosis probability)
Ephyra_Medusa <- function(ii, Sizet1, pars, Esurv, Temp){ 
  # A little fix to prevent the code breaking if and when the model is fed NA temperature entries.
  if(is.na(Temp[ii])){
    return( Esurv[ii] * 0.0095 * 0 * Medusae_size.t(Sizet1[[ii]], pars) )
  } else {
    # this function will use the survival already generated for Ephyra (so that it matches across transitions)
    # probability of producing medusae as a function of Temperature
    No.medusae.mean <- pars["P.medusae.int"] + (Temp[ii] * pars["P.medusae.temp"])
    return( Esurv[ii] * 0.0095 * define_probs(No.medusae.mean, pars["P.medusae.sd"]) * Medusae_size.t(Sizet1[[ii]], pars) )
  }
} # Ishii et al. 2004 showed how the survivorship of Ephyra to metamorphose to Medusae drops to 0.95% of that of survival in newly produced ephyra so that will be replicated here to instigate a bottleneck in Medusae production.

# Define the size of new medusae produced by Ephyra -----------------------------------------------------------
Medusae_size.t <- function(Sizet1, pars){
  threshold <- pars["min.size"]                        
  medusae_size <- dnorm(Sizet1, mean = threshold, sd = 0.1) # it is nessecary to set some arbitrary variation within the size of new medusae.
  return(medusae_size)
}

# Medusae growth (Not being effected by temperature as there isn't enough data for this) -----------------------------------------------------------
g_z1z <- function(Sizet2, Sizet1, pars){
  mean <- pars["Growth.int"] + (pars["Growth.slope"] * Sizet1)
  sd <- sqrt(pi/2) * exp((pars["Growth.sd.int"] + (pars["Growth.sd.slope"] * Sizet1))) # this allows for the variance in new size to change with initial size
  p_den_grow <- dnorm(Sizet2, mean = mean, sd = sd)
  return(p_den_grow)
}

# Convert medusae size into a wet weight estimate (for subsequent estimation of larval production as a function of size) -----------------------------------------------------------
Size_WW_convert <- function(Sizet1, pars){
  # length weight relationship is WW = aL^b
  a <- rnorm(1, pars["a"], pars["a.sd"])
  b <- rnorm(1, pars["b"], pars["b.sd"])
  # run equation
  WW <- a * Sizet1^b
  # return output
  return(WW)
}

# Number of larvae produced as a function of Medusae size -----------------------------------------------------------
# This is being constrained by known physiological thresholds of A. aurita
M_larvae <- function(WWt1, pars, Temp, Sal){
  # define mean larval output as a function of size
  mean <- pars["No.Plan.int"] + (pars["No.Plan.slope"] * WWt1)
  # define variance
  sd <- pars["No.Plan.sd"]
  # estimate output for model
  no_larvae <- exp(rnorm(1, mean, sd))
  # A little fix to ensure reproduction doesn't occur below certain temperatures
  if (Temp < pars["rep.min.temp"] || Sal < pars["rep.min.sal"]) { no_larvae = 0 }
  # and return
  return(no_larvae)
}

# Construct Medusae survival and growth kernel -----------------------------------------------------------
Pz1z <- function (Sizet2, Sizet1, pars, survM) { # This function is fed a predefined survival probability to ensure survival is consistent across associated growth and fecundity estimates
  # A little fix to prevent the code breaking if and when the model is feed missing temperature entries.
  if(is.na(survM)){
    return( 0 * g_z1z(Sizet2, Sizet1, pars) ) 
  } else {
    return( survM * g_z1z(Sizet2, Sizet1, pars) ) 
  }
} 

# Construct Medusae fecundity kernel -----------------------------------------------------------
Fz1z <- function (Sizet1, pars, Ripe_size, Temp, Sal, survM) { # this function will also be fed a predefined ripe size threshold 
  if(is.na(survM)){
    return( 0 )
  } else {
    # are the individuals larger than the reproductive threshold size 
    if(Sizet1 > Ripe_size) {
      # calculate the wet weight of each medusae
      WW <- Size_WW_convert(Sizet1, pars)
      # using their weight then estimate their larval output
      if(is.na(Temp)) {
        return( 0 ) # no reproduction info if there is not temperature.
      } else {
        return( survM * pars["Sex.ratio"] * M_larvae(WW, pars, Temp, Sal) )
      }
    } else {
      return( 0 ) # no reproduction below the sexual maturity threshold.
    }}
}

# 5. Matrix construction functions called within DMat function -----------------------------------------------------------
# Function for building Panel 1 (Discrete matrix element).
mk_mat1 <- function(ii, mat = mat, dims = dims, index, Plan_stasis, Plan_Polyp, Polyp_Polyp, Ephyra_Ephyra, Polyp_bud, Polyp_Ephyra){
  if(ii >= index){
    # build survival matrix
    if(mat == "U") {
      mat_use <- matrix(data = c(Plan_stasis[ii],0,0,0,
                                 0,Polyp_Polyp[ii],0,0,
                                 Plan_Polyp[ii],0,0,0, # for the purposes of this model there is not survival in the roaming polyp stage
                                 0,0,0,Ephyra_Ephyra[ii]), nrow = dims, ncol = dims, byrow = TRUE) 
    }
    # build clonality matrix
    if(mat == "C") {
      mat_use <- matrix(data = c(0,0,0,0,
                                 0,Polyp_bud[ii],0,0,
                                 0,0,0,0,
                                 0,Polyp_Ephyra[ii],0,0), nrow = dims, ncol = dims, byrow = TRUE) 
    }
    # build fecundity matrix
    if(mat == "F") {
      mat_use <- matrix(0, nrow = dims, ncol = dims, byrow = TRUE) # No fecundity prior to Medusae stage
    }
  } else {
    mat_use <- matrix(as.numeric(NA), nrow = dims, ncol = dims, byrow = TRUE)
  }
  
  # extract output
  return(mat_use)
}

# Function for building  Panel 2 (Discrete to continuous) -----------------------------------------------------------
mk_mat2 <- function(ii, mat = mat, m, dims = dims, index, Plan_Polyp_Medusa_store, Ephyra_Medusa_store){
  if(ii >= index){
    # build survival matrix
    if(mat == "U") {
      mat_use <- matrix(data = c(Plan_Polyp_Medusa_store[[ii]], Plan_Polyp_Medusa_store[[ii]], Plan_Polyp_Medusa_store[[ii]], Ephyra_Medusa_store[[ii]]),
                        nrow = m, ncol = dims, byrow = FALSE) 
    }
    # build clonality and fecunditiy matrices
    if(mat %in% c("C", "F")) {
      mat_use <- matrix(data = c(Plan_Polyp_Medusa_store[[ii]], Plan_Polyp_Medusa_store[[ii]], Plan_Polyp_Medusa_store[[ii]], Plan_Polyp_Medusa_store[[ii]]),                        
                        nrow = m, ncol = dims, byrow = FALSE) 
    }
  } else {
    mat_use <- matrix(as.numeric(NA), nrow = m, ncol = dims, byrow = TRUE)
  }
  
  # extract output
  return(mat_use)
}

# Function for constructing Panel 3 (Continuous Medusae survival and growth) -----------------------------------------------------------            
mk_mat3 <- function(ii, mat = mat, index, h, m, Pz1z_store){
  if(ii >= index){
    # build survival matrix
    if(mat == "U") {
      mat_use <- h[ii] * Pz1z_store[[ii]]
    }
    # build clonality and fecunditiy matrices
    if(mat %in% c("C", "F")) {
      mat_use <- h[ii] * 0 * Pz1z_store[[ii]]
    }
  } else {
    mat_use <- matrix(as.numeric(NA), nrow = m, ncol = m, byrow = TRUE)
  }
  
  # extract output
  return(mat_use)
}

# Function for constructing Panel 4 (Continuous to Discrete) -----------------------------------------------------------
mk_mat4 <- function(ii, mat = mat, m, dims = dims, index, h, Fz1z_store){
  if(ii >= index){
    # build survival and clonality matrices
    if(mat %in% c("C", "U")) {
      mat_use <- matrix(data = c((h[ii] * 0 * Fz1z_store[[ii]]), 
                                 (h[ii] * 0 * Fz1z_store[[ii]]),
                                 (h[ii] * 0 * Fz1z_store[[ii]]),
                                 (h[ii] * 0 * Fz1z_store[[ii]])), nrow = dims, ncol = m, byrow = TRUE)
    }
    # build fecundity matrix
    if(mat == "F") {
      mat_use <- matrix(data = c((h[ii] * Fz1z_store[[ii]]), # the medusae only produce larvae not polyps or ephyra
                                 (h[ii] * 0 * Fz1z_store[[ii]]), 
                                 (h[ii] * 0 * Fz1z_store[[ii]]),
                                 (h[ii] * 0 * Fz1z_store[[ii]])), nrow = dims, ncol = m, byrow = TRUE)
    }
  } else {
    mat_use <- matrix(as.numeric(NA), nrow = dims, ncol = m, byrow = TRUE)
  }
  
  # extract output
  return(mat_use)
}

# Function for constructing matrix arrays -----------------------------------------------------------
build_array <- function(ii, mat = mat,
                        U1, U2, U3, U4,
                        C1, C2, C3, C4,
                        F1, F2, F3, F4) {
  # build survival array
  if(mat == "U") {array_use <- cbind(rbind(U1[[ii]],U2[[ii]]), rbind(U4[[ii]],U3[[ii]]))}
  # build clonality array
  if(mat == "C") {array_use <- cbind(rbind(C1[[ii]],C2[[ii]]), rbind(C4[[ii]],C3[[ii]]))}
  # build fecundity array
  if(mat == "F") {array_use <- cbind(rbind(F1[[ii]],F2[[ii]]), rbind(F4[[ii]],F3[[ii]]))}
  # return output
  return(array_use)
}

# And a function to build a seasonal array - combining survival, clonality and reproduction. -----------------------------------------------------------
mk_matB <- function(ii, matU, matC, matF) {
  array_use <- matU[,,ii] + matC[,,ii] + matF[,,ii]
  # return output
  return(array_use)
}

########################### END OF CODE -------------------------------------------
