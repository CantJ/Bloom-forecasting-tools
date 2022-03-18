# This script outlines the function ValidateForecast which can be used to evaluate the accuracy of the GoJelly bloom prediction tool using observer data.

# The process of validating bloom forecasts will overtime help to highlight gaps in our understanding of Jellyfish population dynamics
# and will help to reveal potential seed locations warranting attention as possible sites of polyp populations. 

# Currently, this function is designed to receive observer data in the form of .csv files detailing the timing and location of jellyfish counts made during cruise surveys
# or .json files detailing the GPS coordinates of sightings input into the GoJelly Jelly Spotter app. 

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
# Date last modified: March 2022
# -----------------------------------------------------------------------------

#----------------------------------------------------
# STEP 1: DEFINE PRIMARY VALIDATION FUNCTION
#----------------------------------------------------

## Primary validation function
ValidateForecast <- function(observeData, # Insert observer data
                             # Next the function needs to be told the format the observeData file is in.
                             # IMPORTANT: This function has been designed to work with observer sightings recorded in the GoJelly JellySpotter app. 
                             # This JellySpotter data is stored in .json files within a unique format, that this function is specifically designed to unpack. 
                             # Selecting the 'json' data file will cause the function to assume it is working with JellySpotter data.
                             # Alternatively, to support the more universal application of this function it has also been designed to accommodate abundance data contained within .csv files.
                             # These csv files must contain the variables of Lat (Latitude), Lon (Longitude), Year, Month, and Popsize (Population size: can be counts or density) for each entry.
                             data.type, # = c("json", "csv")
                             forecastData, # Insert simulated data for validating
                             res = 33, # define raster resolution in km (defaults to 33km)
                             time.only = TRUE, # should validation include an evaluation of temporal predictions ignoring spatial patterns? Defaults to TRUE.
                             spatial.only = TRUE, # should validation include an evaluation of spatial predictions ignoring temporal patterns? Defaults to TRUE.
                             # define the spatial details of the forecast data (defaults to the Baltic Region)
                             lonmx = 31, lonmn = 3, latmx = 66, latmn = 53, 
                             # Tell the function the year with which the forecasted data corresponds
                             # IMPORTANT: there must be observer entries also corresponding with this selected year.
                             Year.select){
  
  # Load package dependencies
  packages <- c("rjson", "purrr", "stringr", "raster", "sf", "ggsn", "viridis", "ggplot2", "rnaturalearth",
                "rnaturalearthdata", "scales", "gmodels", "dplyr", "pbapply", "gganimate", "gifski")
  lapply(packages, require, character.only = TRUE)

## 1. Ensure observer data meets desired format requirements 
  
  # The function applies differing approaches if it is dealing with excel data or JellySpotter app data
  if(data.type == "csv"){
    
    # If dealing with csv data this should have already been formatted correctly. Here the function will check if this is the case.
    # Include a little fix to ensure the datafile contains the relevant variables.
    if(any(!is.element(c("Lat", "Lon", "Year", "Month", "Popsize"), names(observeData)))) { stop("Data needs to include each of the following named variables: Lat (Latitude), Lon (Longitude), Year, Month, and Popsize (Population size: can be counts or density)") } 
    
  }
  
  # If dealing with JellySpotter data the function needs to unpack the complex data format to reveal the nessecary details.
  if(data.type == "json"){
  
    # Firstly, clarify temporal and spatial details
    observeData$Month <- observeData$Year <- observeData$Lat <- observeData$Lon <- NA
    # what years are present within the JellySpotter data?
    # the JellySpotter data is recorded in chronological order so the first and last entries will confirm the temporal range of the data
    dims <- dim(observeData)
    YearMin <- as.numeric(paste0(substr(observeData$date[1], start = 8, stop = 12))); YearMax <- as.numeric(paste0(substr(observeData$date[dims[1]], start = 8, stop = 12)))
    # Define an month and year index
    Month_index <- str_sub(month.name, start = 1, end = 3)
    Year_index <- YearMin:YearMax
    # loop through date information extracting required month and year details
    for(ii in 1:length(Month_index)){
      observeData[grep(Month_index[ii], observeData$date), "Month"] <- month.name[ii]
    }
    for(ii in 1:length(Year_index)){
      observeData[grep(Year_index[ii], observeData$date), "Year"] <- Year_index[ii]
    }
  
    # Secondly, unpack GPS location details
    observeData$Lat <- as.numeric(str_sub(observeData$location, 6, regexpr(",", observeData$location)-1))
    observeData$Lon <- as.numeric(str_sub(observeData$location, regexpr(",", observeData$location)+7, nchar(observeData$location)))
  
    # Next, assign a numerical code to quantify the size of observed jellyfish populations 
    # In the JellySpotter app observers can select categories indicating the 'number' of jellyfish observed: Only one < Some < a lot. 
    # Here, these categories are assigned a numerical value to represent this scale 1 = Only one to 3 = a lot. NOTE: 0 will be used to represent no jellyfish with -1 used to indicate no observer effort.
    # This approach is important for maintaining consistency across categorical counts and continuous density estimates. 
    observeData$Popsize <- NA
    Size_index <- c("Only One", "Some", "A Lot")
    for(ii in 1:length(Size_index)){
      observeData[which(observeData$amount == Size_index[ii]), "Popsize"] <- ii
    }
    # To improve the quantity of usable JellySpotter data we can assume that records of none Aurelia individuals correspond with the observation of no Aurelia.
    observeData[which(observeData$chosenJelly != "Moon jelly"),]$Popsize <- 0 # 0 and -1 will be used to differentiate between no individuals observed (0) and no observations carried out (-1)

    # Finally, remove unnecessary variables
    observeData$amount <- NULL
    observeData$chosenJelly <- NULL
    observeData$date <- NULL
    observeData$location <- NULL
    observeData$imageUploadPath <- NULL
    observeData$uuid <- NULL
    
  }
  
  # Regardless of supplied data types - the observeData file will now contain the nessecary information. 
  # Ensure the required variables are in the desired format
  observeData$Lat <- as.numeric(paste(observeData$Lat))
  observeData$Lon <- as.numeric(paste(observeData$Lon))
  observeData$Month <- as.factor(observeData$Month)
  observeData$Year <- as.factor(observeData$Year)
  observeData$Popsize <- as.numeric(paste(observeData$Popsize))
  
  # The data is ready to go!

## 2. Snip observer data to correspond with the temporal and spatial extent of the forecasted data
  
  # remove data from unrequired years
  observeData <- observeData[which(observeData$Year == Year.select),]

  # We are only interested in observations made within the user defined region
  observeData <- observeData[which(observeData$Lat <= latmx & observeData$Lat >= latmn),]
  observeData <- observeData[which(observeData$Lon <= lonmx & observeData$Lon >= lonmn),]

  # Determine the spatial extent of the observer data.
  # This defines a GPS bounding box for the validation
  lonmx2 <- max(observeData$Lon, na.rm = TRUE)
  lonmn2 <- min(observeData$Lon, na.rm = TRUE)
  latmx2 <- max(observeData$Lat, na.rm = TRUE)
  latmn2 <- min(observeData$Lat, na.rm = TRUE)
  
  # Importantly it is nessecary here to determine the range in population size contained within the observer data. 
  # This step is nessecary to scale the forecasted population densities to match observer recording format to accommodate for comparing between forecasted and observed data.
  # For instance the categorical data from the JellySpotter app contains a population size range of 1 - 3 - it is nessecary for the forecasted data to be re-scaled to this spectrum in order to compare the two datasets.
  # This step of the function allows the function to be implemented across categorical and continuous observer data
  # max population size observed?
  maxObs <- max(observeData$Popsize, na.rm = TRUE)
  # smallest non-zero density observed?
  minObs <- min(observeData$Popsize[which(observeData$Popsize != 0)], na.rm = TRUE)
  
  
## 3. Determine spatial and temporal sighting patterns
  
  # Generate a series of spatial rasters showing the distribution of spotter records over time
  
  ## Spatial and temporal patterns combined -------------------------------
  
  # Produce monthly observation rasters
  raster_store <- lapply(1:12, ObserverRaster, observeData = observeData, res = res, xmx = lonmx2, xmn = lonmn2, ymx = latmx2, ymn = latmn2)
  
  ## Spatial patterns only (if requested) ----------------------------------
  if(spatial.only == TRUE){
    
    # determine resolution adjustment
    # there is roughly 111km per degree of latitude and longitude
    raster.dims <- floor(111/res)
    
    # generate blank raster template
    All_r <- raster(nrow = ((latmx2 - latmn2)*raster.dims), ncol = ((lonmx2 - lonmn2)*raster.dims), # this plots at a resolution of ~33km (matching the Risk Map resolution)
                    xmx = lonmx2, xmn = lonmn2, ymx = latmx2, ymn = latmn2)
    # and fill with empty values
    raster::values(All_r) <- NA
    
    # extract spatial details
    xcoords <- matrix(observeData$Lon, ncol = 1)
    ycoords <- matrix(observeData$Lat, ncol = 1)
    # and 'density' data
    PopCount <- matrix(observeData$Popsize, ncol = 1)
    # condense together
    dat <- data.frame(Lon = xcoords, Lat = ycoords, density = PopCount)
  
    # convert to spatial points file
    dat <- st_as_sf(dat, coords = c("Lon", "Lat"), crs = 4326)
    # and generate a density raster file 
    All_r <- rasterize(dat, All_r, 'density', fun = max, background = -1)
    
  }
  
  ## Temporal patterns only (if requested) ----------------------------------
  # This part doesn't produce a raster, instead it visualizes the monthly variation in observed jellyfish densities (spatial variables ignored)
  if(time.only == TRUE){
    
    # Produce plot
    TempPlot <- ggplot() +
      geom_boxplot(aes(x = Month, y = Popsize, fill = Month), data = observeData, width = 0.5, color = "gray", alpha = 0.4) +
      scale_x_discrete(limits = month.name, breaks = month.name, labels = str_sub(month.name, start = 1, end = 3), expand = c(0,0)) +
      scale_fill_viridis(discrete = TRUE, option = "D", direction = 1) +
      scale_y_continuous(expand = c(0, 0), 
                         breaks = c(0, floor(seq(0,maxObs, length.out = 4))[-c(1,4)], maxObs), 
                         labels = c(0, floor(seq(0,maxObs, length.out = 4))[-c(1,4)], maxObs)) + # sets the axis to density appropriate values density
      theme_classic() +
      theme(legend.position="none") +
      theme(legend.title=element_blank()) +
      theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1)) +
      theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.6)) +
      theme(axis.text.y = element_text(size = 12)) +
      xlab("Month") +
      ylab("Medusae Population size") +
      theme(axis.title = element_text(size = 15))
    
  }
  
## 4. Determine forecasted spatial and temporal patterns and convert into raster format
  
  # Separate relevant forecast data components
  Med_den_mean <- forecastData[["Dens"]]
  GPS_lat <- forecastData[["Lat"]]
  GPS_lon <- forecastData[["Lon"]]
  # define length variables
  n_sites <- length(Med_den_mean)
  n_month <- 12
  
  # condense together days from corresponding months - across all requested sites.
  Med_mean <- lapply(1:n_sites, array_reformat, initial_array = Med_den_mean, n_month = n_month)
  Lat <- lapply(1:n_sites, array_reformat, initial_array = GPS_lat, n_month = n_month)
  Lon <- lapply(1:n_sites, array_reformat, initial_array = GPS_lon, n_month = n_month)
  
  # Now ensure that each matrix only contains entries for only times for which there is corresponding GPS coordinates and visa versa.
  Med_mean[is.na(Lat)] <- NA
  Lat[is.na(Med_mean)] <- NA
  Lon[is.na(Med_mean)] <- NA
  
  ## Spatial and temporal patterns combined -------------------------------
  
  # generate rasters
  print("Generating forecast rasters") # progress notifier
  site_store <- pblapply(1:n_sites, ForecastRaster, 
                         Nmax = maxObs, Nmin = minObs, 
                         Lat = Lat, Lon = Lon, Med_mean = Med_mean, 
                         res = res, xmx = lonmx2, xmn = lonmn2, ymx = latmx2, ymn = latmn2, ignore.time = FALSE)
  # separate out raster files and data
  site_raster_store <- lapply(site_store, '[', 'rast')
  if(time.only == TRUE) { site_data_store <- lapply(site_store, '[', 'df') }
  # remove primary storage to save memory
  rm(site_store)
  
  ## Spatial patterns only (if requested) ----------------------------------
  if(spatial.only == TRUE){
    
    # generate rasters
    print("Generating forecast rasters (ignoring temporal patterns)") # progress notifier
    site_store2 <- pblapply(1:n_sites, ForecastRaster, 
                            Nmax = maxObs, Nmin = minObs,
                            Lat = Lat, Lon = Lon, Med_mean = Med_mean, 
                            res = res, xmx = lonmx2, xmn = lonmn2, ymx = latmx2, ymn = latmn2, ignore.time = TRUE)
    
    # Condense rasters across sites into one raster. Importantly here we are only interested in whether the model predicts jellyfish appearance in any given location or not.
    # Therefore if multiple sites predict occurrence for a given location, subsequent rasters can mask previous entries - there is no need to take a total.
    # tell r how to deal with overlaps
    site_store2$fun <- max
    # and condense
    site_r <- do.call(raster::mosaic, site_store2)
    # remove primary storage to save memory
    rm(site_store2)
    
  }
  
  ## Temporal patterns only (if requested) ----------------------------------
  # FIRST MOMENT OF TRUTH
  # DO TEMPORAL EXPECTIONS MATCH (IGNORING SITES)?
  # Here the function will overlay the temporal patterns predicted in the forecast data onto the observed temporal patterns.
  # This will allow for visualizing the correlation between temporal predictions and observations
  if(time.only == TRUE){
    
    # Condense together datafiles across sites
    PredDens <- do.call(rbind, sapply(site_data_store, '[[', 'df'))
    PredDens$Month <- as.factor(PredDens$Month) # generate monthly factor
    
    # extract monthly summaries
    suppressWarnings( TempData <- PredDens %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(mean = mean(DensScale, na.rm = TRUE),
                sd = sd(DensScale, na.rm = TRUE),
                ci.lower = ci(DensScale, na.rm = TRUE)[2],
                ci.upper = ci(DensScale, na.rm = TRUE)[3],
                se = ci(DensScale, na.rm = TRUE)[4]) )
    
    # Overlay temporal trends in predicted densities on observer density plot)
    TempPlot <- TempPlot +
      geom_ribbon(aes(x =  Month, ymin = mean - sd, ymax = mean + sd, group = 1), data = TempData, fill = "blue", colour = NA, alpha = 0.15) +
      geom_line(aes(x = Month, y = mean, group = 1), data = TempData, color = "blue", linetype = "longdash", alpha = 0.5) +
      geom_point(aes(x =  Month, y = mean), color = "blue", data = TempData, size = 2, alpha = 0.5) +
      coord_cartesian(ylim = c(0, maxObs))
    
  }
 
## 5. DETERMINE THE DEGREE OF ERROR ACROSS PREDICTED RASTERS
  # SECOND MOMENT OF TRUTH
  # this will be done by overlaying the observer data rasters with the forecasted rasters and quantifying the margin of error across months.
  # the margin of error will be quantified using Root Mean Square Error (RMSE).
  
  # The function has already produced a plot comparing the degree of overlap in the timing of peak medusae densities observed both in the observer data and the forecasted data (if requested)
  # All that is left to compute is the extent to which the observer and forecasted density rasters agree.
  
  # generate RMSE estimates across months for each site
  print("Validating rasters") # progress notifier
  site_RMSE_store <- pblapply(1:n_sites, RasterValidate, raster_store = raster_store, site_raster_store = site_raster_store)
  
  # extract the outputs for easier visualization/manipulation (here the RMSE estimate will be matched to the GPS coordinates of their corresponding release sites
  rmseData <- as.data.frame(matrix(unlist(site_RMSE_store), ncol = 1, byrow = TRUE))
  # add in nessecary details
  rmseData$Month <- rep(month.name, n_sites); rmseData$Lat <- rep(sapply(Lat, '[[', 1), length.out = dim(rmseData)[1]); rmseData$Lon <- rep(sapply(Lon, '[[', 1), length.out = dim(rmseData)[1])
  names(rmseData) <- c("RMSE", "Month", "Lat", "Lon")
  rmseData$Month <- as.factor(rmseData$Month)
  
  ## Temporal patterns only (if requested) ----------------------------------
  # How does RMSE change across months (The temporal accuracy of predictions)
  if(time.only == TRUE){
  # Determine mean RMSE statistics across the year 
  suppressWarnings( RMSEData1 <- rmseData %>%
                      dplyr::group_by(Month) %>%
                      dplyr::summarise(mean = mean(RMSE, na.rm = TRUE),
                                       sd = sd(RMSE, na.rm = TRUE), 
                                       ci.lower = ci(RMSE, na.rm = TRUE)[2],
                                       ci.upper = ci(RMSE, na.rm = TRUE)[3]) )
  # determine RMSE estimate range
  y.min <- floor(min(RMSEData1$ci.lower, na.rm = TRUE))
  y.max <- ceiling(max(RMSEData1$ci.upper, na.rm = TRUE))
  
  # Produce plot
  plot1 <- ggplot(aes(x = Month, y = mean), data = RMSEData1) +
    geom_ribbon(aes(x = Month, ymin = mean - sd, ymax = mean + sd, group = 1), fill = "blue", alpha = 0.15) +
    geom_line(aes(x = Month, y = mean, group = 1), color = "blue", linetype = "longdash", alpha = 0.5) + 
    scale_y_continuous(limits = c(y.min,y.max)) +
    scale_x_discrete(limits = month.name, breaks = month.name, labels = str_sub(month.name, start = 1, end = 3)) +
    scale_fill_viridis(discrete = TRUE, option = "D", direction = 1) +
    theme_classic() +
    theme(legend.position="none") +
    theme(legend.title=element_blank()) +
    theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1)) +
    theme(axis.text.x = element_text(size = 12, angle = 30, vjust = 0.6)) +
    theme(axis.text.y = element_text(size = 12)) +
    xlab("Month") +
    ylab("Root Mean Square Error (RMSE)") +
    theme(axis.title = element_text(size = 15))
  
  }
  
  ## Spatial and temporal patterns combined --------------------------------
  # How do the monthly accuracy of these predictions varies with site
  
  # format and fix order of month variable 
  rmseData$Month = factor(rmseData$Month, levels = month.name)
  rmseData <- rmseData[order(rmseData$Month),]
  # determine RMSE range
  minRMSE <- min(rmseData$RMSE, na.rm = TRUE)
  maxRMSE <- max(rmseData$RMSE, na.rm = TRUE)
  
  # extract world shapefile
  world <- ne_countries(scale = "medium", returnclass = "sf")
  # define color scale
  pal <- c(viridis(1000, direction = -1, option = "magma"))

  # Build base plot ready for animation
  Base_plot <- ggplot() + 
    geom_sf(data = world, fill = "gray79", colour = "gray47") + 
    geom_point(aes(x = Lon, y = Lat, color = RMSE, group = Month), size = 4, data = rmseData) +
    coord_sf(xlim = c(lonmn, lonmx), ylim = c(latmn, latmx), expand = FALSE) +
    scale_color_gradientn(colours = pal,
                          limits = c(minRMSE, maxRMSE),
                          breaks = c(minRMSE, maxRMSE),
                          labels = c(round(minRMSE,2), round(maxRMSE, 2)), 
                          na.value = NA,
                          guide = guide_colorbar(barwidth = 15, "cm")) + 
    theme_classic() +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) +
    theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1)) +
    theme(legend.text = element_text(size = 20)) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(plot.title = element_text(size = 20)) +
    ggsn::scalebar(data = world, transform = TRUE, dist = 200, dist_unit = "km", model = "WGS84",
                   height = 0.002, location = "bottomright", st.dist = 0.002, st.bottom = FALSE, st.size = 0.1, anchor = c(x = 29.0, y = 53.5))
  
  # Add animation features
  map_with_animation <- Base_plot +
    transition_states(Month) +
    enter_fade() +
    ggtitle('RMSE: {next_state}')
  
  # set up animation
  print("Rendering animated plot")
  suppressWarnings( plot2 <- animate(map_with_animation, fps = 3,
                                     height = 800, width = 900,
                                     renderer = gifski_renderer(loop = F)) )

  ## Spatial patterns only (if requested) ----------------------------------
  if(spatial.only == TRUE){
    
    # Estimate the RMSE between the two single (sites and time combined) rasters.
    AllRMSE <- RMSE(values(All_r), values(site_r))
    
  }
  
  # end of validation function
  # the outputs returned depend on the requested outputs
  if(time.only == TRUE & spatial.only == TRUE) { return(list(RMSE_all_sites = rmseData, # Temporal patterns in prediction error across all sites
                                                             Total_RMSE = AllRMSE, # Prediction error pooling across time and space
                                                             TempPlot1 = TempPlot, # How do density predictions map over temporal patterns in observer sightings
                                                             TempPlot2 = plot1, # How does predictive error vary across time (sites pooled)
                                                             RMSE_map = plot2)) } # animated plot showing how the predictive error of simulations originating across each site vary over time (i.e. which sites produce the best representation of observer sightings)
  if(time.only == TRUE & spatial.only == FALSE) { return(list(RMSE_all_sites = rmseData, Total_RMSE = NA, TempPlot1 = TempPlot, TempPlot2 = plot1, RMSE_map = plot2)) }
  if(time.only == FALSE & spatial.only == TRUE) { return(list(RMSE_all_sites = rmseData, Total_RMSE = AllRMSE, TempPlot1 = NA, TempPlot2 = NA, RMSE_map = plot2)) }
  if(time.only == FALSE & spatial.only == FALSE) { return(list(RMSE_all_sites = rmseData, Total_RMSE = NA, TempPlot1 = NA, TempPlot2 = NA, RMSE_map = plot2)) }
}


#----------------------------------------------------
# STEP 2: DEFINE INTERNAL FUNCTIONS CALLED BY PRIMARY VALIDATION FUNCTION
#----------------------------------------------------

## Define rasterising function
ObserverRaster <- function(ii, 
                          observeData,
                          # user defined resolution
                          res,
                          # spatial bounding box of observer data
                          xmx, xmn, ymx, ymn){
  
  # determine resolution adjustment
  # there is roughly 111km per degree of latitude and longitude
  raster.dims <- floor(111/res)
  
  # extract relevant observer data
  data_use <- observeData[which(observeData$Month == month.name[ii]),]
  
  # if there is no observer data for the selected month - produce a blank raster
  if(dim(data_use)[1] == 0) {
    r <- raster(nrow = ((ymx - ymn)*raster.dims), ncol = ((xmx - xmn)*raster.dims), 
                xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
    # fill it with empty values
    raster::values(r) <- -1  # = no observer effort
  } else {
    # otherwise.... 
    # generate blank raster template
    r <- raster(nrow = ((ymx - ymn)*raster.dims), ncol = ((xmx - xmn)*raster.dims), 
                xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
    # fill it with empty values
    raster::values(r) <- NA
  
    # extract spatial details
    xcoords <- matrix(data_use$Lon, ncol = 1)
    ycoords <- matrix(data_use$Lat, ncol = 1)
    # and 'density' data
    PopCount <- matrix(data_use$Popsize, ncol = 1)
    # condense together
    dat <- data.frame(Lon = xcoords, Lat = ycoords, density = PopCount)
  
    # convert to spatial points file
    dat <- st_as_sf(dat, coords = c("Lon", "Lat"), crs = 4326)
    # and generate a density raster file 
    r <- rasterize(dat, r, 'density', background = -1) # background = no observer effort
  }
    
  # return monthly observation raster
  return(r)
}
  
## Define function for reformatting arrays to combine observations across corresponding days. 
array_reformat <- function(rp, initial_array, n_month){
  tmp_array <- asplit(initial_array[[rp]][,,], MARGIN = 1) # split apart and condense together corresponding matrix rows from across the array
  new_array <- sapply(1:n_month, function(x){ t(tmp_array[[x]])}, simplify = 'array') # transpose aggregated rows
  return( new_array )
}

## Define re-scaling function (to match projected density units with spotter categories)
DensRescale <- function(Oi, Nmax, Nmin, Omin, Omax){
  if(is.na(Oi)) { 
    Rv <- NA 
  } else {
    if(Oi == 0){
      Rv <- 0 # a zero corresponds with no observed individuals and must remain constant
    } else { 
      Rv <- (Nmax - Nmin)/(Omax - Omin)*(Oi - Omax) + Nmax
    }
  }
  return(Rv)
}  
  
## Define function for producing rasters from forecast data
ForecastRaster <- function(ii,  
                           Nmax, Nmin,
                           Lat, Lon, Med_mean,
                           # user defined resolution
                           res,
                           # calculate spatial patterns only? (ignoring temporal aspects).
                           ignore.time = FALSE,
                           # spatial bounding box of observer data
                           xmx, xmn, ymx, ymn){
  
  # determine resolution adjustment
  # there is roughly 111km per degree of latitude and longitude
  raster.dims <- floor(111/res)
  
  # define required parameters
  n_month <- 12
  n_day <- 30
  
  if(ignore.time == FALSE) { # This will work through each month separately for each site 
  
    # generate internal storage
    df <- rast <- list()
  
    # run loop
    for(tt in 1:n_month) {
      # define raster variables
      ycoords <- matrix(Lat[[ii]][,tt,], ncol = 1)
      xcoords <- matrix(Lon[[ii]][,tt,], ncol = 1)
      PopDens <- matrix(Med_mean[[ii]][,tt,]/n_day, ncol = 1)
      # determine scaled density predictions (i.e convert risk into the 1,2,3,4 format used in the spotter records)
      PopDensScaled <- sapply(PopDens, DensRescale, Nmax = Nmax, Nmin = Nmin, Omin = min(PopDens, na.rm = TRUE), Omax = max(PopDens, na.rm = TRUE))
      
      # combine variables together and convert to a spatial points file
      dat <- data.frame(Lon = xcoords, Lat = ycoords, Density = PopDens, DensScale = PopDensScaled)
      dat <- dat[complete.cases(dat),]
      datSP <- st_as_sf(dat, coords = c("Lon", "Lat"), crs = 4326)
      
      # generate a density raster file (Note the raster is being defined for the observer region of interest only)
      r <- raster(nrow = ((ymx - ymn)*raster.dims), ncol = ((xmx - xmn)*raster.dims), # this plots at a resolution of ~33km (same as spotter data rasters)
                  xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
      # and rasterise the scaled density vector
      r <- rasterize(datSP, r, 'DensScale', background = 0) # the difference between observer and forecast data is that in forecast data, no entry means no abundance.
      
      # associate density predictions with corresponding month of the year
      dat$Month <- month.name[tt]
      
      # snip the dataframe to remove entries that do not correspond with the spotter data bounding region
      dat <- dat[which(dat$Lon <= xmx & dat$Lon >= xmn),]
      dat <- dat[which(dat$Lat <= ymx & dat$Lat >= ymn),]
      
      # and return outputs
      df[[tt]] = dat
      rast[[tt]] = r
    }
  
    # return outputs
    return(list(df = df, rast = rast))
  
  } else { # otherwise pool all projected densities across months.
    
    # define raster variables
    ycoords <- matrix(Lat[[ii]], ncol = 1)
    xcoords <- matrix(Lon[[ii]], ncol = 1)
    PopDens <- matrix(Med_mean[[ii]]/n_day, ncol = 1)
    # determine scaled density predictions (i.e convert risk into the 1,2,3,4 format used in the spotter records)
    PopDensScaled <- sapply(PopDens, DensRescale, Nmax = 4, Nmin = 2, Omin = min(PopDens, na.rm = TRUE), Omax = max(PopDens, na.rm = TRUE))
    
    # combine variables together and convert to a spatial points file
    dat <- data.frame(Lon = xcoords, Lat = ycoords, DensScale = PopDensScaled)
    dat <- dat[complete.cases(dat),]
    datSP <- st_as_sf(dat, coords = c("Lon", "Lat"), crs = 4326)
    
    # generate a density raster file (Note the raster is being defined for the observer region of interest only)
    r <- raster(nrow = ((ymx - ymn)*raster.dims), ncol = ((xmx - xmn)*raster.dims), # this plots at a resolution of ~33km (same as spotter data rasters)
                xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
    # and rasterise the scaled density vector
    r <- rasterize(datSP, r, 'DensScale', fun = max, background = 0) # background for the forecasted data = no predicted density.

    # return outputs
    return(rast = r)
  }
}  

## Define a function that will estimate the root mean squared error for paired rasters.
RMSE <- function(x, y) { 
  
  # Firstly in the current set up comparing observer rasters with no data (all 0) with forecast rasters in distorted and biased RMSE estimates.
  # Therefore it is not wise to compare forecasted rasters when observer rasters consist of no observations.
  if( all(x == -1) ) {
    RMSE1 <- NA
  } else {
    # Estimate RMSE for only rasters in which observer raster contain observer effort.
    RMSE1 <- sqrt(mean((x - y)^2, na.rm = TRUE))
  
    ## quick tidy up
    if(is.nan(RMSE1)) { RMSE1 <- NA }
  }
  
  # and return outputs
  return(RMSE1 = RMSE1)
} 

## Define function for applying validation across site rasters
RasterValidate <- function(ii, raster_store, site_raster_store){
  
  # raster_store contains the observer data (each raster corresponding with monthly observations)
  # split out monthly observations
  Jan <- raster_store[[1]]; Feb <- raster_store[[2]]; Mar <- raster_store[[3]]; Apr <- raster_store[[4]]; May <- raster_store[[5]]; Jun <- raster_store[[6]]
  Jul <- raster_store[[7]]; Aug <- raster_store[[8]]; Sep <- raster_store[[9]]; Oct <- raster_store[[10]]; Nov <- raster_store[[11]]; Dec <- raster_store[[12]]
  
  # now lay out predicted monthly rasters
  Jan2 <- site_raster_store[[ii]]$rast[[1]]; Feb2 <- site_raster_store[[ii]]$rast[[2]]; Mar2 <- site_raster_store[[ii]]$rast[[3]]; Apr2 <- site_raster_store[[ii]]$rast[[4]] 
  May2 <- site_raster_store[[ii]]$rast[[5]]; Jun2 <- site_raster_store[[ii]]$rast[[6]]; Jul2 <- site_raster_store[[ii]]$rast[[7]]; Aug2 <- site_raster_store[[ii]]$rast[[8]]
  Sep2 <- site_raster_store[[ii]]$rast[[9]]; Oct2 <- site_raster_store[[ii]]$rast[[10]]; Nov2 <- site_raster_store[[ii]]$rast[[11]]; Dec2 <- site_raster_store[[ii]]$rast[[12]]
  
  # and determine the margin of error between corresponding monthly rasters
  # if the RMSE is returned as NaN due to the abundance of NAs across some rasters, then the function will instead return NA to make the outputs easier to deal with.
  em_Jan <- RMSE(values(Jan), values(Jan2))
  em_Feb <- RMSE(values(Feb), values(Feb2))
  em_Mar <- RMSE(values(Mar), values(Mar2))
  em_Apr <- RMSE(values(Apr), values(Apr2))
  em_May <- RMSE(values(May), values(May2))
  em_Jun <- RMSE(values(Jun), values(Jun2))
  em_Jul <- RMSE(values(Jul), values(Jul2))
  em_Aug <- RMSE(values(Aug), values(Aug2))
  em_Sep <- RMSE(values(Sep), values(Sep2))
  em_Oct <- RMSE(values(Oct), values(Oct2))
  em_Nov <- RMSE(values(Nov), values(Nov2))
  em_Dec <- RMSE(values(Dec), values(Dec2))
  
  # and return outputs
  return(list(em_Jan, em_Feb, em_Mar, em_Apr, em_May, em_Jun, em_Jul, em_Aug, em_Sep, em_Oct, em_Nov, em_Dec))
}

########################### END OF CODE -------------------------------------------
