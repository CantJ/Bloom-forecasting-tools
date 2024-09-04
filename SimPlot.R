# This script outlines the SimPlot function used to plot forecasted densities provided to it as part of the GoJelly predictive tool for forecasting the occurrence of jellyfish blooms.

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
# -----------------------------------------------------------------------------

#----------------------------------------------------
# STEP 1: DEFINE FUNCTION
#----------------------------------------------------

# Plotting function that will produce visual outputs using details returned by the JellySim function. 
SimPlot <- function(meanRast, confRast, sites, xmx, xmn, ymx, ymn) {
  
  # Load package dependencies
  packages <- c("terra", "sf", "viridis", "ggplot2", "rnaturalearth",
                "rnaturalearthdata", "scales", 'plyr')
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    invisible(install.packages(packages[!installed_packages]))
  }
  invisible(lapply(packages, library, character.only = TRUE))

  # Estimate maximum recorded densities
  dmax <- c(round_any(max(values(meanRast)), 0.01, f = ceiling),
            round_any(max(values(confRast)), 0.01, f = ceiling))
  
  # read in world shape file
  world <- ne_countries(scale = "medium", returnclass = "sf") 
  
  # For plotting, the raster layers will need coercing into stacked data-frames
  # Extract coordinates
  plotCoords <- xyFromCell(meanRast, seq_len(ncell(meanRast)))
  # Extract density estimates
  meanValues <- stack(as.data.frame(values(meanRast)))
  confValues <- stack(as.data.frame(values(confRast)))
  # reassign variable names
  names(meanValues) <- names(confValues) <- c('Density', 'Month')
  levels(meanValues$Month) <- levels(confValues$Month) <- paste0('Month ', 1:nlyr(meanRast))
  colnames(plotCoords) <- c('Lon', 'Lat')
  # And combine coordinates and values back together
  meanDat <- cbind(plotCoords, meanValues)
  confDat <- cbind(plotCoords, confValues)
  
  # format release site gps coordinates for plotting
  sites <- st_as_sf(sites, coords = c("Lon", "Lat"), crs = 4326) # this crs code corresponds with a WGS84 projection
  
  # define color scale
  pal <- c("white", viridis(1000, direction = -1, option = "magma"))
  # Define facet labels
  monthNames <- as_labeller(c("Month 1" = month.abb[1], "Month 2" = month.abb[2], "Month 3" = month.abb[3], "Month 4" = month.abb[4],
                            "Month 5" = month.abb[5], "Month 6" = month.abb[6], "Month 7" = month.abb[7], "Month 8" = month.abb[8], 
                            "Month 9" = month.abb[9], "Month 10"= month.abb[10], "Month 11" = month.abb[11], "Month 12" = month.abb[12], "Month 13" = month.abb[1]))
  
  # Build plots
  # Mean density forecasts
  meanMap <- ggplot() +
    geom_tile(aes(x = Lon, y = Lat, fill = Density), data = meanDat) +
    facet_wrap(~Month, labeller = monthNames) +
    geom_sf(data = world, fill = "gray79", colour = "gray47") + 
    geom_sf(data = sites, fill = "black", size = 3.5, shape = 21) + # add location of selected release sites to the plot
    coord_sf(xlim = c(xmn, xmx), ylim = c(ymn, ymx), expand = FALSE) +
    scale_fill_gradientn(colours = pal,
                         breaks = c(0, dmax[1]),
                         labels = c('Low', 'High'),
                         limits = c(0, dmax[1]), 
                         guide = guide_colorbar(barwidth = 15, "cm", ticks = F)) + 
    # define axis tick frequency
    scale_y_continuous(breaks = c(ymn+1,  ymx-2)) +
    scale_x_continuous(breaks = c(xmn+2,  xmx-2)) +
    theme_classic() +
    theme(legend.position="top", legend.title=element_blank(), legend.background = element_rect(fill = '#EFF0F8')) +
    theme(axis.line.x = element_line(linewidth = 1), axis.line.y = element_line(linewidth = 1)) +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
    theme(legend.text = element_text(size = 15)) +
    theme(strip.text = element_text(size = 15, colour = '#EFF0F8')) +
    theme(panel.spacing = unit(1, "lines")) +
    theme(plot.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8', linewidth = 0), 
          panel.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8')) +
    theme(strip.background = element_rect(fill = '#2C2D7C')) +
    xlab("\nLongitude") +
    ylab("Latitude\n") +
    theme(plot.title = element_text(size = 20))
  
  # Forecast confidence
  confMap <- ggplot() +
    geom_tile(aes(x = Lon, y = Lat, fill = Density), data = confDat) +
    facet_wrap(~Month) +
    geom_sf(data = world, fill = "gray79", colour = "gray47") + 
    geom_sf(data = sites, fill = "black", size = 3.5, shape = 21) +
    coord_sf(xlim = c(xmn, xmx), ylim = c(ymn, ymx), expand = FALSE) +
    scale_fill_gradientn(colours = pal,
                         breaks = c(0, dmax[2]),
                         labels = c('Low', 'High'),
                         limits = c(0, dmax[2]), 
                         guide = guide_colorbar(barwidth = 15, "cm", ticks = F)) + 
    # define axis tick frequency
    scale_y_continuous(breaks = c(ymn+1,  ymx-2)) +
    scale_x_continuous(breaks = c(xmn+2,  xmx-2)) +
    theme_classic() +
    theme(legend.position="top", legend.title=element_blank(), legend.background = element_rect(fill = '#EFF0F8')) +
    theme(axis.line.x = element_line(linewidth =  1), axis.line.y = element_line(linewidth = 1)) +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
    theme(legend.text = element_text(size = 15)) +
    theme(strip.text = element_text(size = 15, colour = '#EFF0F8')) +
    theme(panel.spacing = unit(1, "lines")) +
    theme(plot.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8', linewidth = 0),
          panel.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8')) +
    theme(strip.background = element_rect(fill = '#2C2D7C')) +
    xlab("\nLongitude") +
    ylab("Latitude\n") +
    theme(plot.title = element_text(size = 20))
  
  # return maps
  return(list(mean = meanMap, confidence = confMap))
}

########################### END OF CODE -------------------------------------------
