# This script outlines the SimPlot function used to plot forecasted densities provided to it as part of the GoJelly predictive tool for forecasting the occurrence of jellyfish blooms.

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
# Date last modified: March 2022
# -----------------------------------------------------------------------------

#----------------------------------------------------
# STEP 1: DEFINE FUNCTION
#----------------------------------------------------

# Plotting function that will produce visual outputs using details returned by the JellySim function. 
SimPlot <- function(mean_plots, conf_plots, dmax, confmax, rel_coords, xmx, xmn, ymx, ymn) {
  
  # Load package dependencies
  packages <- c("raster", "sf", "ggsn", "viridis", "ggplot2", "rnaturalearth",
                "rnaturalearthdata", "scales")
  lapply(packages, require, character.only = TRUE)
  
  # read in world shape file
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  # take the plot stores and generate raster stacks consisting of all monthly rasters
  mean_plot_stack <- raster::stack(mean_plots)
  conf_plot_stack <- raster::stack(conf_plots)
  
  # For plotting, the raster stacks will need coercing into stacked data-frames
  # Extract coordinates
  mean_coords <- xyFromCell(mean_plot_stack, seq_len(ncell(mean_plot_stack)))
  conf_coords <- xyFromCell(conf_plot_stack, seq_len(ncell(conf_plot_stack)))
  # Extract density estimates
  mean_values <- stack(as.data.frame(getValues(mean_plot_stack)))
  conf_values <- stack(as.data.frame(getValues(conf_plot_stack)))
  # reassign variable names
  names(mean_values) <- names(conf_values) <- c('value', 'Month')
  levels(mean_values$Month) <- levels(conf_values$Month) <- month.name
  colnames(mean_coords) <- colnames(conf_coords) <- c('long', 'lat')
  # And combine coordinates and values back together
  mean_data <- cbind(mean_coords, mean_values)
  conf_data <- cbind(conf_coords, conf_values)
  
  # format release site gps coordinates for plotting
  rel_coords <- st_as_sf(rel_coords, coords = c("Lon", "Lat"), crs = 4326) # this crs code corresponds with a WGS84 projection
  
  # define color scale
  pal <- c("white", viridis(1000, direction = -1, option = "rocket"))
  
  # Build plots
  # Mean density forecasts
  mean_map <- ggplot() +
    geom_tile(aes(x = long, y = lat, fill = value), data = mean_data) +
    facet_wrap(~Month) +
    geom_sf(data = world, fill = "gray79", colour = "gray47") + 
    geom_sf(data = rel_coords, fill = "black", size = 3.5, shape = 21) + # add location of selected release sites to the plot
    coord_sf(xlim = c(xmn, xmx), ylim = c(ymn, ymx), expand = FALSE) +
    scale_fill_gradientn(colours = pal,
                         breaks = c(0, dmax),
                         labels = c("Low", "High"),
                         limits = c(0, dmax), 
                         guide = guide_colorbar(barwidth = 15, "cm")) + 
    # define axis tick frequency
    scale_y_continuous(breaks = c(ymn+1,  ymx-2)) +
    scale_x_continuous(breaks = c(xmn+2,  xmx-2)) +
    theme_classic() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.background = element_rect(fill = '#EFF0F8')) +
    theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1)) +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    theme(legend.text = element_text(size = 15)) +
    theme(strip.text = element_text(size = 15, colour = '#EFF0F8')) +
    theme(panel.spacing = unit(1, "lines")) +
    theme(plot.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8', size = 0), 
          panel.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8')) +
    theme(strip.background = element_rect(fill = '#2C2D7C')) +
    xlab("\nLongitude") +
    ylab("Latitude") +
    theme(plot.title = element_text(size = 20)) +
    ggsn::scalebar(data = mean_data, dist = 300, transform = TRUE, dist_unit = "km", model = "WGS84", height = 0.05, location = "bottomright",
                   st.dist = 0.05, st.bottom = FALSE, st.size = 5, anchor = c(x = 29.0, y = 53.5),
                   facet.var = "Month", facet.lev = "December")
  
  # Forecast confidence
  conf_map <- ggplot() +
    geom_tile(aes(x = long, y = lat, fill = value), data = conf_data) +
    facet_wrap(~Month) +
    geom_sf(data = world, fill = "gray79", colour = "gray47") + 
    geom_sf(data = rel_coords, fill = "black", size = 3.5, shape = 21) +
    coord_sf(xlim = c(xmn, xmx), ylim = c(ymn, ymx), expand = FALSE) +
    scale_fill_gradientn(colours = pal,
                         breaks = c(0, confmax),
                         labels = c(0, format(confmax, scientific = FALSE)),
                         limits = c(0, confmax), 
                         guide = guide_colorbar(barwidth = 15, "cm")) + 
    # define axis tick frequency
    scale_y_continuous(breaks = c(ymn+1,  ymx-2)) +
    scale_x_continuous(breaks = c(xmn+2,  xmx-2)) +
    theme_classic() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.background = element_rect(fill = '#EFF0F8')) +
    theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1)) +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
    #theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    theme(legend.text = element_text(size = 15)) +
    theme(strip.text = element_text(size = 15, colour = '#EFF0F8')) +
    theme(panel.spacing = unit(1, "lines")) +
    theme(plot.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8'),
          panel.background = element_rect(fill = '#EFF0F8', color = '#EFF0F8')) +
    theme(strip.background = element_rect(fill = '#2C2D7C')) +
    xlab("\nLongitude") +
    ylab("Latitude") +
    theme(plot.title = element_text(size = 20)) +
    ggsn::scalebar(data = conf_data, dist = 300, transform = TRUE, dist_unit = "km", model = "WGS84", height = 0.05, location = "bottomright",
                   st.dist = 0.05, st.bottom = FALSE, st.size = 5, anchor = c(x = 29.0, y = 53.5),
                   facet.var = "Month", facet.lev = "December")
  
  # return maps
  return(list(mean = mean_map, conf = conf_map))
}

########################### END OF CODE -------------------------------------------
