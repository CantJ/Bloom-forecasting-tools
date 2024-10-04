# GoJelly: Tools to predict and evaluate spatial and temporal patterns in Jellyfish bloom occurrence.
---

This repository is a companion to the manuscript *'Coupling hydrodynamic drifting simulations and seasonal demographics to unmask the drivers of jellyfish blooms'* (in Prep) Cant et al., and the GoJelly Risk Map application (**https://gojellyeu.shinyapps.io/gojellyapp/**). Together the code contained within this repository can be used to predict and visualize spatial and temporal patterns in the aggregation and movements of jellyfish populations.

Our motivation behind the development of this predictive framework was to introduce tools that offer a basis for (1) predicting the appearance of jellyfish blooms and (2), through the comparison of these predictions with recorded observations, resolving gaps in our understanding of jellyfish dynamics. Observational field data on the spatial movements of jellyfish populations is lacking, in part due to their complex and cryptic life-cycle. Ultimately, this limits our understanding of the causes of bloom formation, and makes it difficult to mitigate the economic consequences of jellyfish blooms worldwide. We hope that by making these tools publicly available this work will encourage the development of bloom forecasting, and help to guide the design of monitoring surveys to improve understanding of regional jellyfish populations.

For further details on the various function contained within this repository please check out the corresponding manuscript above, and also don't hesitate to get in touch with **James Cant at james.cant91@gmail.com** with any queries.

This work was carried out as part of the EU Horizon funded project GoJelly. Further information on this project can be found at **https://gojelly.eu/**.



## File Details:

***JellySim.R:***
This script contains the primary forecasting function (and its internal dependent functions) of our bloom forecasting framework. In the manuscript associated with this code repository we have used this function to predicts temporal and spatial patterns in the aggregation of *Aurelia aurita* populations within the Baltic Sea. However, this tool is intentionally designed to be sufficiently flexible for adapting to accommodate differing focal populations and locations. Briefly, the **JellySim** function takes information outlining the location of interest and the hydrodynamic & abiotic conditions associated with that region of interest, and combines these details with a demographic model describing the dynamics of the focal population (see **DMat.R**), to forecast the location, timing, and size, of jellyfish aggregations. In essence, the tool can be provided alternative details regarding location, abiotic conditions, and demographic details, to modify it taxonomic and spatial focus.

***DMat.R:***
This script contains all of the demographic functions needed to parameterise a hybrid matrix population model simultaneously describing the discrete and continuous dynamics of *A. aurita* populations within the Baltic Sea. In the context of our forecasting framework, this script serves as a demographic black box, and deals with all the population modelling elements of our predictive tool. This approach ensures that it is possible for users to redefine and parameterise their own demographic functions specific to their study organism, which can then be used to modify the scope of the forecasting functions in **JellySim.R**. ***Note***: if modifying *DMat* the details contained with **Demographic parameters.csv** need to be modified accordingly.

***SimPlot.R:***
This script deals with plotting the outputs returned by **JellySim.R**. Provided the spatial details provided to **SimPlot.R** correspond with those provided to **JellySim.R** this plotting function can be adapted to any desired location.

***ParSens.R***
In lieu of sufficient data for carrying out formal model validation this script contains a function for testing the sensitivity of modelled outputs to the specific input parameters. We anticipate that this sensitivity tool will facilitate the identification of key underlying processes within the modeled dynamics of jellyfish populations that warrant additional attention during future surveys.

***app.R:***
This script contains the UI and server functions underpinning the interactive Risk Map application currently being hosted on shinyapps.io (*https://gojellyeu.shinyapps.io/gojellyapp/*). Ultimately, this script combines the simulating, modelling, and visualisation functions, contained with the **JellySim.R**, **PeriodicMat.R**, and **SimPlot.R** scripts into a more user-friendly format.

***ParticleTracks2021.RData:***
This data file contains the hydrodynamic (drifting pathways) and abiotic conditions (Temperature & Salinity) experienced by drifting individuals within the Baltic Sea region, following there release at various locations and times throughout 2021.

***Release_site_GPS_info:***
This csv file contains GPS coordinates corresponding with the release sites used in producing **ParticleTracks2021.R**. It is not required in the context of the forecasting tools, but is used by **app.R** to generate visuals.

***Demographic parameters.csv:***
A csv file containing modelling coefficients used to parameterise the vital rate functions contained within **DMat.R**. These specific coefficients have been estimated using empirical data obtained from *A. aurita* populations within the Baltic region.
