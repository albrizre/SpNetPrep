library(shiny)
library(shinythemes)
library(leaflet)
# library(colourpicker)
# library(leaflet.minicharts)
# library(shinyWidgets)
# library(rhandsontable)
# library(ggmap)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("yeti"),
                  
  navbarPage(HTML(paste(strong("SpNetPrep: Facilitating spatial statistics on road networks"))), selected="Network Edition", id="nav", 
            
    tabPanel("Network Edition",
             
      h1(HTML(paste("Spatial Linear Network",strong("Edition")))),
             
        sidebarLayout(
          sidebarPanel(width = 2,
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                        h4("INSERT A ROAD NETWORK")
            ),           
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.5em",
              fluidRow(
                fileInput('file', 'Upload a SpatialLines or SpatialLinesDataFrame (.RDS)',
                          accept=c('.Rdata'))
                )
              ),
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-2.5em;text-align:center",
                h4("____________________________")
            ), 
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                h4("EDIT YOUR NETWORK")
            ),    
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em",
                fluidRow(
                  column(6, sliderInput('vertex_radius','Vertex Radius',min=1,max=10,value=5,step=1)), 
                  column(6, sliderInput('edge_thickness','Edge Thickness',min=1,max=10,value=5,step=1))
                )
            ),
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-1.5em;margin-top:0.25em;text-align:left",
              fluidRow(
              radioButtons(inputId="action", label="Choose an action:",
                            choiceNames = list("Join vertex","Remove edge","Add point (+edge)","Add two points (+edge)"),
                            choiceValues =  list("add.edge","remove.edge","add.point","add.twopoints"))
              )
            ),
                            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-1.5em;text-align:center",
                h4("____________________________")
            ), 
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                h4("SIMPLIFY YOUR NETWORK")
            ),
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0.2em",
                fluidRow(
                  column(6, textInput(inputId="angle_max",label="Angle:",
                            value="", width = "100%")),
                  column(6,textInput(inputId="length_max",label="Length:",
                            value="", width = "100%"))
                )
            ),
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-bottom:1.5em",
                fluidRow(
                  actionButton("simplifynet", "Simplify linear network", class = "btn-primary", style="width: 100%;")
                )
            ),
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-1em;text-align:center",
                h4("____________________________")
            ), 
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em; margin-top:-0.8em;text-align:center",
                h4("EXECUTE AND EXPORT")
            ),  
            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0.5em;",
              fluidRow(
                actionButton("rebuildnet", "Rebuild linear network", class = "btn-primary", width="100%")
              )
            ),
                            
            div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em",
              fluidRow(
                downloadButton("downloadnet", "Download linear network", class = "btn-primary", style="width: 100%;")
              )
            )
          ),
               
        mainPanel(
          leafletOutput("mymap",height=730,width=1345)
        )
      )
    ),
   
   tabPanel("Network Direction",
            
            h1(HTML(paste("Spatial Linear Network",strong("Direction")))),
            
            sidebarLayout(
              sidebarPanel(width = 2,
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                               h4("INSERT A ROAD NETWORK")
                           ),           
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.5em",
                               fluidRow(
                                 fileInput('fileNetworkDirection', 'Upload a SpatialLines or SpatialLinesDataFrame (.RDS)',
                                           accept=c('.Rdata'))
                               )
                           ),
                           
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-2.5em;text-align:center",
                               h4("____________________________")
                           ), 
                           
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                               h4("ADD DIRECTION TO YOUR NETWORK")
                           ),                
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-1.5em",
                               fluidRow(
                                  radioButtons(inputId="action_flow", label="Choose an action:",
                                              choiceNames = list("Add flow","Add long flow","Remove flow","Remove long flow"),
                                              choiceValues =  list("add.flow","add.longflow","remove.flow","remove.longflow"))
                               )
                           ),
                           
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-1em;text-align:center",
                               h4("____________________________")
                           ), 
                           
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em; margin-top:-0.8em;text-align:center",
                               h4("EXECUTE AND EXPORT")
                           ),  
                           
                           
                           div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em",
                               fluidRow(
                                 downloadButton("downloadnetNetworkDirection", "Download linear network", class = "btn-primary", style="width: 100%;")
                               )
                           )
              ),
              
              mainPanel(
                  leafletOutput("mymapdirection",height=730,width=1345)
              )
            )
   ),
   
  tabPanel("Point Pattern Revision",
             
             h1(HTML(paste("Spatial Point Pattern on a Linear Network",strong("Revision")))),
             
             sidebarLayout(
               sidebarPanel(width = 2,
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                      h4("Insert a point pattern on a road network")
                  ),  
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:1.5em",
                      fluidRow(
                        selectInput("utm_select", "Select UTM zone:", c(1:60), selected = NULL, multiple = FALSE,
                                    selectize = TRUE, width = NULL, size = NULL)
                      )
                  ),
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0.5em; margin-top:-2.5em",
                      fluidRow(
                        fileInput('filePointPattern_lppx', 'Upload a lpp + ppx object (.RDS format)',
                                  accept=c('.Rdata'))
                      )
                  ),
                
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-2.5em;text-align:center",
                      h4("____________________________")
                  ), 
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-0.2em; margin-top:-0.8em;text-align:center",
                      h4("VISUALIZE AND REVISE YOUR NETWORK")
                  ),  
                            
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em",
                    fluidRow(
                       radioButtons("clusterevents_edition", "Cluster events",
                                    choices = c(Yes = "Yes",
                                                No = "No"),
                                                selected = "No")
                    )
                  ),
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:-1.5em",
                      fluidRow(
                        radioButtons(inputId="action_point", label="Choose an action:",
                                     choiceNames = list("Explore pattern","Move event"),
                                     choiceValues =  list("explore.pattern","move.event"))
                      )
                  ),
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:2em; margin-top:-1em;text-align:center",
                      h4("____________________________")
                  ), 
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em; margin-top:-0.8em;text-align:center",
                      h4("EXECUTE AND EXPORT")
                  ),  
                  
                  
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0.5em;",
                      fluidRow(
                        actionButton("rebuildpoint", "Rebuild point pattern", class = "btn-primary", width="100%")
                      )
                  ),
                  div(style = "font-size: 15px; padding: 0px 0px; margin-bottom:0em",
                      fluidRow(
                        downloadButton("downloadpointpattern", "Download point pattern", 
                                       class = "btn-primary", style="width: 100%;")
                      )
                    )
                  ),
               mainPanel(
                 leafletOutput("mymapPointPatternEdition",height=730,width=1345)
               )
             )
    ),
    
    tabPanel("ReadMe", withMathJax(), includeMarkdown("ReadMe.Rmd"))
  )
))
