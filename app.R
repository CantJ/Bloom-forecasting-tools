# This script underpins the GoJelly interactive Jellyfish bloom and aggregation forecasting tool.
# Firstly the script defines all nessecary modelling functions (also contained separately within the JellySim.R, DMat.R and SimPlot.R script files)
# then the script outlines the ui and server functions needed for compiling the web app.

# Primary Author: James Cant
# Contact: james.cant91@gmail.com
# -----------------------------------------------------------------------------

#----------------------------------------------------
# STEP 1: DEFINE PACKAGE AND DATA REQUIREMENTS
#----------------------------------------------------

# Load required packages
require('spsUtil', quietly = TRUE)
require('data.table', quietly = TRUE)
require('terra', quietly = TRUE)
require('sf', quietly = TRUE)
require('purrr', quietly = TRUE)
require('pbapply', quietly = TRUE)
require('abind', quietly = TRUE)
require('foreach', quietly = TRUE)
require('future', quietly = TRUE)
require('doFuture', quietly = TRUE)
require('purrr', quietly = TRUE)
require('gamlss', quietly = TRUE)
require('gamlss.dist', quietly = TRUE)
require("viridis", quietly = TRUE)
require("ggplot2", quietly = TRUE)
require("rnaturalearth", quietly = TRUE)
require("rnaturalearthdata", quietly = TRUE)
require("scales", quietly = TRUE)
require('plyr', quietly = TRUE)
require('ggsn', quietly = TRUE)
require('shiny', quietly = TRUE)
require('shinydashboard', quietly = TRUE)
require('shinyWidgets', quietly = TRUE)
require('shinythemes', quietly = TRUE)
require('plotly', quietly = TRUE)

# load required data
load("ParticleTracks2021.RData") # Hydrodynamic information
params <- read.csv("Demographic parameters.csv", row.names = 1) # Vital rate coefficients
params <- with(params, setNames(x, rownames(params)))
sites <- read.csv("Release_site_GPS_info.csv") # release location coordinates (corresponds with drifting pathways)

# define model parameters that are not user defined.
# simulation iteration interval
nd = 30
# simulation spatial extents (lat/long)
xx = 31
xn = 3
yx = 66
yn = 53
# Demographic matrix dimensions
m = 200
# Multicore processing control
pl = FALSE # parallel processing not possible on shiny server

#----------------------------------------------------
# STEP 2: LOAD SIMULATION MODEL FUNCTION DEPENDANCIES
#----------------------------------------------------

source("DMat.R") # demographic functions
source("JellySim.R") # forecasting tools
source("SimPlot.R") # plotting functions

#----------------------------------------------------
# STEP 3: GENERATE INTERACTIVE RELEASE POINT MAP
#----------------------------------------------------

# convert release location GPS datafile into a shape file for plotting
siteCoords <- st_as_sf(sites, coords = c("Lon", "Lat"), crs = 4326) # this crs code corresponds with a WGS84 projection
# and for plotting purposes re add the latitude and longitude information to this data file
siteCoords$Lat <- sites$Lat; siteCoords$Lon <- sites$Lon
# read in world shape file
world <- ne_countries(scale = "medium", returnclass = "sf") 

# and plot the release sites.
sitePlot <- suppressWarnings(ggplot() +
  geom_sf(data = world, fill = "gray79", colour = "gray47") + 
  geom_sf(aes(fill = Region, col = Region, names = Index, Lat = Lat, Lon = Lon), data = siteCoords, size = 1.8, shape = 21) +
  coord_sf(xlim = c(xn, xx), ylim = c(yn, yx), expand = FALSE) +
  scale_fill_manual(values = c("black","navy","sandybrown","pink","lightskyblue1","dodgerblue","burlywood3", "cornsilk2"),
                    guide = guide_legend(override.aes = list(size = 10, shape = 22))) +
  scale_color_manual(values = c("black","navy","sandybrown","pink","lightskyblue1","dodgerblue","burlywood3", "cornsilk2")) +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1)) +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 15)) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggsn::scalebar(data = world, transform = TRUE, dist = 200, dist_unit = "km", model = "WGS84",
                 height = 0.004, location = "bottomright", st.dist = 0.002, st.bottom = FALSE, anchor = c(x = 32.5, y = 53.5)))

#----------------------------------------------------
# STEP 4: DEFINE WEB APP SPECIALIST FUNCTIONS
#----------------------------------------------------

# define User interface 
ui <- fluidPage(
  navbarPage(
    id = "PageNavigation",
    # Main title
    title = "GoJelly Bloom Risk Map",
    collapsible = TRUE, 
    inverse = FALSE, # light text on a dark background?
    theme = shinytheme("superhero"),
    # Firstly insert a home page introducing users to the application.
    tabPanel(strong("Home"),
             setBackgroundColor(
               color = "#EFF0F8",
               gradient = "linear",
               direction = "bottom",
               shinydashboard = FALSE),
             # Insert some visual details that will remain consistent across the application
             tags$head(
               tags$style(
                HTML('
                  .navbar { box-shadow: 0 0 15px rgba(0, 0, 0, 0.5); }
                  .navbar .navbar-nav {float: left; 
                        color: #EFF0F8; 
                        font-size: 15px;
                        background-color: #2C2D7C ;}
                  .navbar.navbar-default.navbar-static-top{ color: #EFF0F8; 
                                            background-color: #2C2D7C ;}
                  .navbar-default .navbar-nav > .active > a, 
                      .navbar-default .navbar-nav > .active > a:focus, 
                      .navbar-default .navbar-nav > .active > a:hover {
                                      color: #2C2D7C;
                                      background-color: #DC9257;} 
                  .dropdown-menu { color: #EFF0F8; 
                                   background-color: #2C2D7C ;}             
                  .dropdown-menu > .active > a {color: #2C2D7C;
                                                background-color: #DC9257;}
                  .dropdown-menu > li > a:hover{color: #EFF0F8;
                                      background-color: #DC9257;}
                  .nav .open > a, .nav .open > a:focus, .nav .open > a:hover {
                                                                background-color: #DC9257 !important;}
                  .tabbable > .nav > li > a {
                      color:  #2C2D7C ;
                      background-color: #EFF0F8 ; }
                  .tabbable > .nav > li a:hover, .tabbable > .nav > .active > a {
                      color: #DC9257 !important;
                      background-color: #2C2D7C !important; }
                  .navbar-default .navbar-nav > li > a:hover{ color: #EFF0F8;
                                      background-color: #DC9257;} 
                  .navbar .navbar-header {float: left; } 
                  .navbar-default .navbar-brand { color: #EFF0F8; 
                                            width: 500px;
                                            height: 50px;
                                            font-size: 35px;
                                            text-align: center;
                                            background-color: #2C2D7C ;}
                  #sidebar { background-color: #2C2D7C; }
                  h3{
                    font-family: "Arial";
                    color: #EFF0F8;
                  }
                  h6 {
                    font-family: "Arial";
                    color: #EFF0F8;
                  }
                  h4 {
                    color: #ACAFD9;
                    font-family: "Arial";
                  }
                  h5 {
                    color: #2C2D7C;
                    font-family: "Arial";
                  }
                  .modal-footer{
                  background-color: #000000;
                  }
                  .modal-body{
                  background-color: #000000;
                  }
                  .footer {
                   position: auto;
                   bottom: 0;
                   width: 100%;
                   height: 200px;
                   background-color: #000000; 
                   box-shadow: 0 0 15px rgba(0, 0, 0, 0.5);
                  }')),
              tags$script(" $(document).ready(function () {
                  $('#start').on('click', function (e) {
                  window.scrollTo(0, 0)
                    });
                    });")
                  ),
             fluidRow(
               column(9,
                 h3(strong("Welcome"), style = 'color: #2C2D7C; text-align: center;'),
                 h5(p("Jellyfish are an important component of marine ecosystems worldwide, contributing to the transfer of nutrients between coastal and oceanic environments and therefore supporting a diverse range of food webs and biodiversity. 
                      Many jellyfish species even present novel sources for various medicinal products.
                      Yet, alongside their ecological value, the sudden appearance of vast numbers of jellyfish during so-called bloom events can also be detrimental to the fishing, aquaculture, and tourism industries."),
                    br(),
                    p("In recent decades the occurrence of problematic jellyfish blooms is thought to have increased as a consequence of ecosystem degradation.
                      However, it is still debated whether this increase in human-jellyfish conflict is a result of growing numbers of jellyfish, or simply because we are continually expanding how we are using the worlds oceans."),
                    br(),
                    p("Evaluating patterns in the timing and location of jellyfish blooms is therefore crucial for not only mitigating their economic impacts, but also for understanding the status of global jellyfish populations.")),
                  style = 'text-align: center;'),
                 br(),
                 br(),
                column(2, a(img(src = "GoJelly_logo.png", style = "height: 250px; display: block; margin-left: auto; margin-right: auto;"), href = "https://gojelly.eu/", target = "_blank",
                            br()))),
             fluidRow(img(src = 'Photo_banner.png', width = '1400px', style = "display: block; margin-left: auto; margin-right: auto; border-radius: 10px; box-shadow: 0 0 15px rgba(0, 0, 0, 0.3) "))
             ),
    tabPanel(strong("Objectives"),
             fluidRow(column(2), column(8, h3(p(strong("Why is the Risk Map important?")), style = 'color: #2C2D7C; text-align: center')), column(2)),
             fluidRow(column(4), 
                      column(4, 
                             h5(p(em("'All models are wrong, but some are useful'.")),
                                style = 'text-align: center;'),
                             h5(p(strong("Geroge Box"),"(1987)", em("British Statistician.")),
                                style  = 'text-align: right')), 
                      column(4)),
             fluidRow(
               column(2),
               column(8, h5(p("We cannot ignore the impact uncertainty has on our ability to predict ecological events, but it should never prevent us from making predictions.
                              Making predictions, and then validating them, offers the opportunity to identify what it is about ecological systems we do not yet understand, and therefore provides insight into where further work is required.
                              Using predictive tools to first highlight and then resolve gaps in our knowledge is just as valuable as their use in anticipating ecological events."),
                            br(),
                            p("By linking the spatial movements, demographic profiles, and potential release locations, of jellyfish populations this Risk Map tool will help to enhance understanding of jellyfish populations and bloom formation:"), style  = 'text-align: center'),
                         h5(p(HTML(paste0("<ul><li>", strong(em("As the location(s) of polyp beds are revealed, this Risk Map will support exploration into where subsequent blooms can be expected")), "</li>"))),
                            p(HTML(paste0("<li>", strong(em("This tool will also enable the matching of observed blooms with potential seed locations revealing as of yet undocumented polyp populations.")), "</li></ul>"))),
                            br(),
                         style  = 'text-align: left')),
               column(2)),
             fluidRow(img(src = 'Photo_banner.png', width = '1400px', style = "display: block; margin-left: auto; margin-right: auto; border-radius: 10px; box-shadow: 0 0 15px rgba(0, 0, 0, 0.3)"))
             ),
    # Now a tab to give information on Jellyfish, GoJelly and the interactive map
    navbarMenu(strong("About"),
               tabPanel("European jellyfish",
                        fluidRow(
                          column(5, h3(strong("Common European jellyfish"), style = 'color: #2C2D7C; text-align: center'),
                        h5(p(HTML(paste0("Jellyfish is a collective term used to refer to over 1000 species of gelatinous (‘jelly-like’) marine organisms. 
                            In Europe the most common jellyfish species encountered include the Moon jelly (", strong("A"), ", ", em("Aurelia aurita"), ") the Blue jellyfish (", strong("B"), ", ", em("Cyanea lamarckii"), "), the Upside-down jellyfish (", strong("C"), ", ", em("Cassiopea andromeda"), "), the Lion’s mane (", strong("D"), ", ", em("Cyanea capillata"), ") and the Barrel jellyfish (", strong("E"), ", ", em("Rhizostoma pulmo"), "). 
                            The prominence of these species is helped by their recognizable ", em("Scphozoan"), " (true jellyfish) form consisting of a bell-shaped ‘body’ and a series of slender tentacles."))),
                           style = 'text-align: center;'),
                        br(),
                        h5(p(HTML(paste0("Yet, European waters are also home to numerous less-familiar jellyfish species such as the Sea Walnut (", strong("F"), ", ", em("Mnemiopsis leidyi"), "), which is member of the ", em("Ctenophora"), " (Comb jellies)."))),
                           br(),
                           p("Regardless of their differing appearances, however, the majority of jellyfish species possess a capacity for bloom formation. 
                             This bloom formation can occur either directly, following increased reproduction and growth (true blooms), or indirectly, driven by the localised accumulation of multiple populations (apparent blooms)."), 
                          br(),
                          p("Unfortunately, despite the ubiquity of jellyfish species worldwide, the movements and ecology of their populations remain poorly understood. Subsequently, we have limited knowledge for the processes that drive and trigger jellyfish blooms."),
                          style = 'text-align: center')
                        ),
                        column(7, img(src = 'Species_pics.png', width = '800px', style = "display: block; margin-left: auto; margin-right: auto;"),
                                 h5(p(em(strong("Figure 1."), "Photographs of six jellyfish species commonly sighted within European waters.")),
                                    style = "text-align: center;")
                        ))
                       ),
               tabPanel("Jellyfish ecology",
                        fluidRow(
                          column(7, h3(strong("The jellyfish life-cycle"), style = 'color: #2C2D7C; text-align: center'),
                        h5(p("Our limited knowledge for the movements and status of jellyfish populations stems from their complex life-cycle which consists of pelagic (‘drifting’) and benthic (‘fixed to the seabed’) phases (Fig. 2); 
                              a complex life-cycle that is further complicated by high stage-specific plasticity in response to changing environmental conditions – meaning that as jellyfish transition through the cycle their tolerances for different environmental conditions change."),
                           style  = 'text-align: center'),
                        br(),
                        h5(p(strong("1)"), "   Although the number of life stages can differ between species, jellyfish typically start out life as sexually produced", em("planula larvae"), "."),
                           br(),
                           p(strong("2)"), "   After a short period, these free-swimming larvae attach themselves to suitable areas of seabed and develop into polyps, known as", em("scyphistoma"), "."),
                           br(),
                           p(strong("3)"), "   Although, the polyp phase of most jellyfish species is highly cryptic, polyp populations can reach considerable densities as individuals reproduce through asexual reproduction.
                            Equally, it is not completely understood how long different jellyfish species spend as scyphistoma, but it is thought that this phase represents the longest part of the jellyfish life-cycle and can last several years."),
                           br(),
                           p(strong("4)"), HTML(paste0("   Eventually, however, when conditions allow, jellyfish polyps undergo a process known as ", em("strobilation"), ", during which they release free-swimming ", em("ephyra"), "."))),
                           br(),
                           p(strong("5)"), HTML(paste0("   These ephyrae then move away from the substrate, and over time, develop into sexually mature adult ", em("medusae"), ", the characteristic nomadic phase of the jellyfish life-cycle and ultimately the key bloom forming stage."))),
                           style  = 'text-align: left')
                        ),
                        column(5, img(src = 'Life_cycle.png', width = '550px', style = "display: block; margin-left: auto; margin-right: auto;"),
                                 h5(p(em(strong("Figure 2."), "Jellyfish life-cycles are complex and consist of both drifting and sessile benthic phases.")),
                                    style = "text-align: center;")))
                        ),
               tabPanel("The Risk Map",
                        fluidRow(
                          column(5, h3(strong("Technical details"), style = 'color: #2C2D7C; text-align: center'),
                          h5(p("Predicting the timing and location of jellyfish blooms requires following the development of individual jellyfish from their polyp phase into adult medusae.
                               Medusae development and, crucially, the distance they travel, depends on the environmental conditions individuals experience as they drift through the water column.
                               Oceanic currents therefore play a key role in regulating spatial patterns in the movements of jellyfish populations.
                               To recreate this scenario, this Risk Map pairs two components to predict spatial and temporal patterns in the dynamics of jellyfish populations (Fig. 3):"),
                             p(strong("1)"), " A drifting model describing ocean current pathways and the environmental conditions (e.g., temperature and salinity) experienced by drifting particles."), 
                             p(strong("2)"), " A demographic model describing how environmental conditions influence the survival, development, and reproduction of individual jellyfish."),
                             br(),
                             p(HTML(paste0("This Risk Map is currently configured to forecast bloom occurrence in Moon jelly (", em("A. aurita"), ") populations within the Baltic Sea, with efforts to expand these capabilities ongoing.
                               The Risk Map, therefore, contains a series of equations outlining how the prfiles of survival (S) and growth/development (G) vary across the life-cycle of ", em("A. aurita"), ", and how these are influenced by temperature and salinity.
                               The Risk Map then contains a series of pathways followed by drifting particles released within the Baltic Sea, and the environmental conditions they experienced as a result."))),
                             p(HTML(paste0("Selecting from preconfigured release sites, users of the Risk Map can then simulate the release of ephyrae at various coastal locations around the Baltic Sea.
                               The Risk Map then forecasts the development and expected movement of these released individuals, returning heat maps displaying how the density of ", em("A. aurita"), " medusae is expected to vary in both time and space across the Baltic Sea."))),
                             br(),
                             style  = 'text-align: center')),
                          column(7, br(), br(), br(), img(src = 'Work_flow.png', width = '700px', style = "display: block; margin-left: auto; margin-right: auto;"),
                                 br(),
                                 h5(p(em(strong("Figure 3."), "Combining seasonal demographic models describing population specific vital rates with drifting simulations is nessecary to forecast spatial and temporal patterns in the occurrence of jellyfish blooms.")),
                                    style = "text-align: center;"),
                                 br())
                          )
                        ),
               tabPanel("GoJelly",
                        fluidRow(column(4, br(), img(src = "Chrysaora_hysoscella.webp", style = 'width: 400px; display: block; margin-left: auto; margin-right: auto;'),
                                        br()),
                                 column(8, h3(strong("GoJelly: A gelatinous solution to plastic pollution"), style = 'color: #2C2D7C; text-align: center'),
                                        h5(p("GoJelly is a collaborative effort bringing together a range of expertise in technological development, business analysis, fishing practises, and both environmental and social research, from across Denmark, Germany, Norway, Portugal, Slovenia, Israel, Italy, France, and China."),
                                           p("Overall, GoJelly seeks to develop and promote innovative solutions to global marine pollution, using jellyfish as a microplastics filter and as a novel source of blue-green biotechnology products. 
                                              Throughout Europe, jellyfish, particularly when they aggregate in large numbers ('jellyfish blooms'), are typically seen as a nuisance as they can disrupt the fishing and tourism industries.
                                              However, jellyfish can also play a crucial role in our ability to regulate and remove micro- and nano- plastic particles from the world’s oceans. 
                                              When disturbed, many jellyfish species release a mucus that, alongside a variety of functions, serves to remove particulate matter from the water that would otherwise interfere with the jellyfish’s activities. 
                                              It has been shown that this mucus acts as a flocculent and can efficiently concentrate together small plastic particles and even spilled oil. 
                                              Thus, jellyfish mucus represents a potential tool that could enable the collection and removal, and perhaps even the recycling, of plastic pollutants from coastal seas and wastewater effluents."), 
                                           p("The GoJelly consortium is also working to promote jellyfish as a valuable resource for the food and feed, and pharmaceutical industries. 
                                              With this work GoJelly aims to encourage the use of jellyfish as a biological fertilizer for organic farming, whilst also developing jellyfish-based products for use in the development of cosmetics, nutraceutics and medicine."),
                                           style = 'text-align: center'),
                                        br())))
               ),
    # Next insert a tab detailing the user interface components.
    tabPanel(strong("Model Setup"),
             # Insert a sidebar panel to display user defined parameters 
             sidebarLayout(
               sidebarPanel(id = "sidebar",
                            style = "border-radius: 10px; box-shadow: 0 0 15px rgba(0, 0, 0, 0.3)",
                            # Select release location
                            numericInput("rel_location", h5(strong("Select release location(s):"), style = "color: #DC9257"), 564, min = 1, max = 1196),
                            h5(p(style = "text-align: justify; color: #EFF0F8", "Select the location from which to initiate the release of drifting jellyfish individuals. Release locations can be selected by selecting their corresponding index value. Move cursor over the map on the right to display the index values of potential polyp population locations within the Baltic Sea.")),
                            h5(p(style = "text-align: center; color: #DC9257", strong("Note:"), "only one release location can be selected per simulation.")),
                            br(),
                            # define acceptance of default simulation parameters
                            h5(strong("Additional parameters"), style = 'color: #DC9257; text-align: left'),
                            h5(p(style = "text-align: justify; color: #EFF0F8", "In order to forecast spatial and temporal patterns in the aggregation of Jellyfish populations this Risk Map tool requires details on the timing of ephyra release, and the number of times to repeat each simulation in order to account for natural variability in the processes of jellyfish survival and growth.")),
                            br(),
                            h5(em("Ephyra release period:"), style = "color: #DC9257"),
                            fluidRow(
                              column(8, align="center", # Define ephyra release period
                                     sliderTextInput("rel_months", label = NULL,
                                                     choices = c("January", "February", "March", "April", "May", "June", "July",
                                                                 "August", "September", "October", "November", "December"),
                                                     selected = c("February","May"), width = '100%'), 
                                     offset = 2)),
                            br(),
                            h5(em("Number of Iterations:"), style = "color: #DC9257"),
                            fluidRow(
                              column(12, align="center",
                                     # Define number of stochastic iterations 
                                     knobInput("zmax", label = NULL,
                                               min = 2, max = 1000, value = 10,
                                               step = 1,
                                               displayPrevious = TRUE,
                                               fgColor = "#DC9257", 
                                               inputColor = "#DC9257"))),
                            fluidRow(
                              column(10,
                            h5(p(style = "text-align: center", strong("Note:"), "Increasing the number of iterations will extend the time needed for the forecast to run to completion."), style = "color: #DC9257"),
                            br(),
                            offset = 1)),
                            h5(em("See further details tab for more information on each parameter"), style = "color: #EFF0F8"),
                            br(),
                            # Include a start button    
                            actionButton("start", "Run Simulation",
                                         style="color: #EFF0F8; background-color: #C75272; border-color: #fff")
                            ),
               # Next to these user inputs display an interactive map allowing users to selected desired release locations, and a tab detailing the biological relevance of each parameter.
               mainPanel(
                 tabsetPanel(
                   tabPanel("Release locations",
                            h5(p("Move cursor over desired release location to display its Index number and GPS coordinates")),
                            # Insert map
                            box(ggplotly(sitePlot, tooltip = c("names", "Lat", "Lon"), width = 900, height = 750) %>% layout(legend = list(orientation = "h", x = 0, y = -0.2)) %>% config(displayModeBar = F),
                                style = 'display: block; margin-left: auto; margin-right: auto;')
                   ),
                   tabPanel("Further details",
                            wellPanel(
                              h4("Jellyfish release: Locations & timing?", style = 'color: #DC9257; text-align: center'),
                              h5(p("Although people frequently observe drifting medusae, we know very little about where these jellyfish originate. To improve our ability to forecast the movements of jellyfish populations this Risk Map allows exploration into how the existence of polyp populations across various coastal locations could be expected to influence patterns in the appearance of jellyfish medusae. By changing the ",
                                   span(strong("release location"), style = 'color: #DC9257'),
                                   " you can investigate how changing the location from where jellyfish are initially released as ephyra (see life-cycle tab) impacts upon where their populations subsequently aggregate."),
                                 p("Meanwhile, the timing of ephyrae release (strobilation) is regulated by local environmental conditions, and is usually initiated by a sudden drop in temperature followed by seasonal warming. In temperate European seas, strobilation therefore usually occurs between February and May. However, the exact timing of strobilation is highly variable across years and regions, and will be impacted by future climatic shifts in the timing of seasonal temperature cycles. By changing the ",
                                   span(strong("release period"), style = 'color: #DC9257'),
                                   "you can evaluate how changing the timing of ephyra release will influence to subsequent development of medusae populations."),
                                 style = "color: #EFF0F8; text-align: center"),
                              br(),
                              h4("Simulation details", style = 'color: #DC9257; text-align: center'),
                              h5(p("Selecting the",
                                   span(strong("Number of iterations"),style = 'color: #DC9257'),
                                   "informs the Risk Map tool of how many times it needs to repeat simulations. Each time a simulation is initiated, new estimates describing the survival and growth characteristics of the population are calculated.",
                                   "These simulations are, therefore, stochastic simulations, with the Risk Map returning density estimates averaged from across each repeated simulation, plus a measure of the variability across these different simulations (standard deviation). The more iterations (repeats) the Risk Maps completes, the more reliable its projected outputs."),
                                 style = "color: #EFF0F8; text-align: center"), 
                              # change well panel background color
                              style = 'background: #2C2D7C')
                   )))
             )),
    # Next insert a tab to display simulation outputs 
    tabPanel(strong("Results"),
             value = "ResultsPanel",
             # This tab will consist of two sub-tabs - one showing the primary density plots and the other the variance plots
             tabsetPanel(
               # Primary Medusae density plots
               tabPanel("Density patterns",
                        h5(p("This tab displays expected temporal and spatial variation in", 
                             HTML(paste0(strong("relative Medusae density"), ".", ""))),
                           p("Relative Medusae density is displayed using a colour scale depicting the liklihood or",
                             strong("RISK"),
                             "of encounter."),
                           p(em("Note: The solid circles denote the position of the selected release location")),
                           style = "text-align: center"),
                        plotOutput(outputId = 'mean.plot', height = '1300px')
                       ),
               # confidence plot - showing the extent of variation generated by demographic stochasticity.
               tabPanel("Stochastic variability",
                        h5(p(style = "text-align: center", "Plots displayed on this tab show the standard deviation of all projected densities; offering a measure of confidence in any predicted density patterns.")),
                        plotOutput(outputId = 'conf.plot', height = '1300px')
                       )
             )
    ),
    # Insert universal page footer across tabs
    footer = tags$footer(class = "footer",
                         style = "background-color: #000000",
                         p(
                           hr(),
                           column(1, offset = 1, a(img(src = "Logo_dark.png", style = "height: 150px"), href = "https://gojelly.eu/", target = "_blank")),
                           column(3, offset = 1, p(h4(strong("GoJelly"))),
                                  p(h6("A gelatinous solution to plastic pollution: An innovative EU H2020 funded project seeking to use jellyfish as a solution to combat marine litter.")),
                                  p(h6(span(strong("Find out more:"), style = 'color: #ACAFD9'), a("https://gojelly.eu/", href = "https://gojelly.eu/", target = "_blank")))),
                           column(2, p(h4(strong("Get in touch today")), style = "text-align: center"),
                                  p(h6(strong("Jamileh Javidpour:"), "jamileh@biology.sdu.dk"), style = "text-align: center"),
                                  p(h6(strong("Owen Jones:"), "jones@biology.sdu.dk"), style = "text-align: center"),
                                  p(h6(strong("James Cant:"), "jic2@st-andrews.ac.uk"), style = "text-align: center")), 
                           column(3, p(h4(strong("Horizon 2020"))),
                                  p(h6("GoJelly was funded by the European Union’s Horizon 2020 research and innovation programme under grant agreement No 774499.")),
                                  br(),
                                  div("(c) 2018 | GoJelly.eu", style = 'font-size: 12px; text-align: right')))
    )
  )
)

# define server functions (instructions - for how to convert inputs into the outputs
server <- function(input, output, session) { # inputs are the user define variables & outputs are the visuals returned to the user.
  ## create outputs and store - these outputs must use the requested inputs
  
  # Use observeEvent to ensure simulation guides users to outputs once simulation is initiated.
  observeEvent(input$start, {
    updateTabsetPanel(session = session, inputId = "PageNavigation", selected = "ResultsPanel")
  })
  
  # use eventReactive to delay the reactive simulation process until users click 'Run Simulation'.
  # initiating the simulation will call the JellySim function
  sim_data <- eventReactive(input$start, {
    
    # Display loading message (particularly useful for longer simulations)
    showModal(modalDialog(h3("Jellyfish drifting..."),
                          h5(br(),
                             "Hang in there! It can take up to fifteen minutes for a simulation to complete.", style = "color: #EFF0F8; text-align: center"),
                          p(h5("Once this message disappears, please wait a few seconds for the page to update...", style = "color: #EFF0F8; text-align: center")),
                          footer = img(src = "Logo_dark.png", style = "height: 150px"), 
                          size = 'l',
                          easyClose = FALSE,
                          fade = TRUE))
    
    # Takes the categorical month range provided and converts   
    # into the numerical format required by the simulation function
    rel_month_num <- c(match(input$rel_months[1], month.name):match(input$rel_months[2], month.name)) 
    
    # instigate simulation
    Jelly_data <- quiet(JellySim(pars = params, driftData = SINMOD_data, 
                           xmx = xx, xmn = xn, ymx = yx, ymn = yn, # internal data requirements of the function
                           n_days = nd,
                           m = m,
                           # user defined parameters                          
                           rel_location = input$rel_location,
                           rel_months = rel_month_num, 
                           zmax = input$zmax,
                           parallel = pl))
    
    # Update user 
    removeModal()
    
    # let the reactive function know what information it is focusing on
    return(Jelly_data)
  })
  
  # following the implementation of the simulation the SimPlot function will use the required simulation outputs to produce the monthly plots
  sim_outputs <- reactive({
    
    # first determine the required data has been produced.
    req(sim_data())
    
    # create maps
    Jelly_plots <- SimPlot(meanRast= sim_data()$mean, sdRast = sim_data()$conf,
                            sites = sim_data()$site,
                            xmx = xx, xmn = xn, ymx = yx, ymn = yn)
    
    # Now extract each of the required monthly plots
    # Primary Plot
    MeanPlot <- Jelly_plots$mean
    # Confidence Plot
    ConfPlot <- Jelly_plots$conf
    
    # Condense together the plots 
    plots <- list(MeanPlot = MeanPlot, ConfPlot = ConfPlot)
    
    # Let the reactive function know what information it is focusing on
    return(plots)
  })
  
  # Insert desired plot outputs
  # Primary plot
  output$mean.plot <- renderPlot({
    tryCatch( {
      plots <- sim_outputs()
      plots$MeanPlot
    } , shiny.silent.error = function(e) {
      validate("Select release location(s) to initiate density simulation")
    } )
  }, height = 'auto', bg = '#EFF0F8', execOnResize = TRUE)
  
  # Confidence plots
  output$conf.plot <- renderPlot({ 
    tryCatch( {
      plots <- sim_outputs()
      plots$ConfPlot
    } , shiny.silent.error = function(e) {
      validate("")
    } )
  }, height = 'auto', bg = '#EFF0F8', execOnResize = TRUE)
  
  # close server function
} 

# Knit the two elements together
shinyApp(ui = ui, server = server)

## ------------------------------------------------------------------- End of code -----------------------------------------------------------------
