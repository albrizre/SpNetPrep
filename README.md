

# Summary 



The **SpNetPrep** package does not deal with statistics, but with the previous steps that can be required in order to perform a spatial statistical analysis of a point pattern that lies on a linear network representing a road structure. In this regard, the name chosen for the package summarizes its main goal of "Spatial Network Preprocessing" (**SpNetPrep**). The main feature provided by the **SpNetPrep** package is an interactive application that allows carrying out the complete preprocessing of a linear network that comes from a road structure. First, the user needs to install the package via CRAN or via GitHub. Then, the execution of the function *runAppSpNetPrep()* in the R console launches the application allowing its full use, which is also possible to be done online following the link https://albriz.shinyapps.io/spnetprep/. If the application is run from the R console, it is necessary to click the option "Open in browser" when it shows, or define "Run external" for the opening of Shiny applications in order to be able to download the modifications performed on the objects uploaded to it. <br />

According to the technical difficulties that the development of a spatial analysis over a linear network implies, the **SpNetPrep** takes advantage of the R packages **leaflet** and **shiny** to provide an intuitive application that helps to reduce such difficulties. Specifically, the **SpNetPrep** package focuses on the following parts of the preprocessing process that could be required prior to any spatial analyses over a linear network: network creation and edition, network direction endowment and point pattern revision and modification. <br />

In order to start, users can obtain a road network of their interest via the OpenStreetMaps (OSM) platform or from other public or private sources. Hence, when the user is in possession of a road network in a right R format, the **SpNetPrep** application includes a "Network Edition" section that usually would constitute the starting point of the preprocessing phase. At this part of the application, users can introduce their networks in order to delete edges, join vertex to form new edges and create new points that are connected to the preexisting vertex or directly between them. Of course, users that had previously created their road networks with **SpNetPrep** can use this edition section to make changes on them. <br />

The manual edition (or curation) of a linear network representing a road structure is an important step that must be taken in order to correct possible mistakes (not updated road configurations), remove some undesired parts (pedestrian or secondary roads, depending on the application) and also to simplify some zones of the network whose complexity could obscure the analysis being performed (which is sometimes very notorious in round-abouts or complex intersections). Another important question to take into consideration when working with a linear network structure is its directionality. Depending on the kind of dataset being treated, network direction could be of no interest, but this should not be the case when analyzing traffic-related data. For example, in order to use a spatial model with a collection of accident counts at the road segment level (for instance, with the **spdep** package), the provision of a directionality to the linear network would become essential to define a realistic neighbourhood structure. In a similar way, if a geostatistical approach is established (see the **gstat** package) in order to predict a quantitative measure along a road network which is likely affected by traffic flow, the lack of consideration of the directions that can be taken by the vehicles that use it could lead to meaningless results. In this regard, the "Network Direction" section of the **SpNetPrep** application attempts to facilitate the enhancement of a network with this valuable information. <br />

Once the network structure is properly curated and endowed with a direction (if necessary), a point pattern can be located along it from a dataset containing geocoded information. This step can be achieved straight by using the (shortest) orthogonal projection of each pair of coordinates into the linear network, for example with the *project2segment* function of the R package **spatstat**, but depending on the level of accuracy of the coordinates available, it can lead to some percentage of events wrongly placed. Therefore, the "Edit a Point Pattern" section of the application allows the user to investigate this issue while providing a whole picture of the distribution of the point pattern in the road network being studied.

