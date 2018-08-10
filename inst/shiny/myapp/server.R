library(shiny)
library(shinythemes)
library(leaflet)
library(leaflet.minicharts)
library(htmltools)
library(sp)
library(spatstat)
library(maptools)
library(geosphere)
library(prodlim)
library(rgeos)
library(raster)
library(rgdal)

### AUXILIARY FUNCTIONS

VerticesFormanEje <- function(grafo, eje){
  tabla=cbind(grafo$vertices$x,grafo$vertices$y)
  inicio=c(grafo$lines$ends$x0[eje],grafo$lines$ends$y0[eje])
  fin=c(grafo$lines$ends$x1[eje],grafo$lines$ends$y1[eje])
  buscarinicio=row.match(inicio,tabla)
  buscarfin=row.match(fin,tabla)
  buscar=c(buscarinicio,buscarfin)
  return(buscar)
}

VerticesEjes <- function(grafo){
  vertices_ejes=c()
  for (i in c(1:grafo$lines$n)){
    print(i)
    vertices_ejes=rbind(vertices_ejes,VerticesFormanEje(grafo,i))
  }
  return(vertices_ejes)
}

AnguloPuntos <- function(grafo,origeni,fini){
  resta=fini-origeni
  anguloi=atan(resta[2]/resta[1])*(180/pi)
  return(anguloi)
}

BuscarVertice <- function(V,Punto,d){
  dist_x=abs(V$x-Punto$x)
  dist_y=abs(V$y-Punto$y)
  buscar=which((dist_x<d & dist_y<d)==T)
  return(buscar)
}

DetectarEjeGrafo <-function(grafo,datos){
  ejes=c()
  for (j in c(1:length(datos$x))){
    #print(j)
    x=datos$x[j]
    y=datos$y[j]
    aux=list(x=x,y=y)
    X_aux=lpp(aux,grafo)
    eje=X_aux$data$seg
    ejes=c(ejes,eje)
  }
  return(ejes)
}

UTM2LONLAT <- function(coord,zone){
  if (class(coord)=="numeric"){
    coord=t(as.matrix(coord))
  }
  df=as.data.frame(coord)
  colnames(df)=c("lon","lat")
  coordinates(df) <- c("lon", "lat")
  proj4string(df) <- CRS(paste0("+proj=utm +zone=",zone," ellps=WGS84"))
  res <- spTransform(df, CRS("+proj=longlat +datum=WGS84"))
  res=as.data.frame(cbind(res@coords[,1],res@coords[,2]))
  colnames(res)=c("lon","lat")
  return(res)
}

LONLAT2UTM <- function(coord,zone){
  if (class(coord)=="numeric"){
    coord=t(as.matrix(coord))
  }
  df=as.data.frame(coord)
  colnames(df)=c("lon","lat")
  coordinates(df) <- c("lon", "lat")
  proj4string(df) <- CRS("+proj=longlat +datum=WGS84")
  res <- spTransform(df, CRS(paste0("+proj=utm +zone=",zone," ellps=WGS84")))
  res=as.data.frame(cbind(res@coords[,1],res@coords[,2]))
  colnames(res)=c("lon","lat")
  return(res)
}

shortestpath <- function(L, i, j) {
  L <- as.linnet(L, sparse=FALSE)
  d <- L$dpath
  m <- L$m
  to <- L$to
  from <- L$from
  path <- i
  leader <- i
  repeat {
    k <- setdiff(which(m[leader,]), path)
    leader <- k[which.min(d[i,k] + d[k, j])]
    path <- c(path, leader)
    if(leader == j) break
  }
  return(path)
}

LongFlowBetweenPoints <- function(grafo, v_origin, v_end, Vert, Vertices){
  ok=F
  iter=0
  while (ok==F & iter<=50){
    ### Identify vertex connected with origin
    search=which((Vert[,1]==v_origin | Vert[,2]==v_origin)==T)
    connected_vertex=unique(as.numeric(Vert[search,]))
    connected_vertex=connected_vertex[connected_vertex!=v_origin]
    print(connected_vertex)
    flow=c()
    for (i in c(1:length(connected_vertex))){
      angle=AnguloPuntos(grafo,c(Vertices$x[v_origin],Vertices$y[v_origin]),c(Vertices$x[connected_vertex[i]],Vertices$y[connected_vertex[i]]))
      if (abs(angle)<=5 | abs(180-angle)<=5){
        flow=c(flow,connected_vertex[i])
        if (connected_vertex[i]==v_end){
          ok=T
        } else{
          iter=iter+1
          print(iter)
        }
      }
    }
  }
  return(flow)
}

SimplificarRed <- function(grafo,angulo_max,longitud_max,vertices_ejes_aux){
  grado_2=which((vertexdegree(grafo)==2)==T)
  for (i in c(1:length(grado_2))){
    #print(i)
    ### Se toma el correspondiente v?rtice de grado 2 de los originales
    ### Se hallan los dos ejes que contienen a este v?rtice
    ### Se hallan los dos puntos que son unidos con v mediante los ejes
    v=grado_2[i]
    buscar_1=which((vertices_ejes_aux[,1]==v)==T)
    buscar_2=which((vertices_ejes_aux[,2]==v)==T)
    ejes_eliminar=c(buscar_1,buscar_2)
    puntos_unir=c(vertices_ejes_aux[buscar_1,],vertices_ejes_aux[buscar_2,])
    puntos_unir=puntos_unir[puntos_unir!=v]
    
    ### Obtener longitudes y ?ngulos ejes
    
    ### ?ngulos ejes 
    
    resta1=c(grafo$vertices$x[puntos_unir[1]],grafo$vertices$y[puntos_unir[1]])-
      c(grafo$vertices$x[v],grafo$vertices$y[v])
    angulo1=atan(resta1[2]/resta1[1])*(180/pi)
    
    resta2=c(grafo$vertices$x[puntos_unir[2]],grafo$vertices$y[puntos_unir[2]])-
      c(grafo$vertices$x[v],grafo$vertices$y[v])
    angulo2=atan(resta2[2]/resta2[1])*(180/pi)
    
    ### Longitudes ejes 
    
    longitud1=pointDistance(c(grafo$vertices$x[v],grafo$vertices$y[v]),
                            c(grafo$vertices$x[puntos_unir[1]],grafo$vertices$y[puntos_unir[1]]),
                            lonlat = F)
    longitud2=pointDistance(c(grafo$vertices$x[v],grafo$vertices$y[v]),
                            c(grafo$vertices$x[puntos_unir[2]],grafo$vertices$y[puntos_unir[2]]),
                            lonlat = F)
    
    condicion_angulo=abs(angulo1-angulo2)<=angulo_max | (180-abs(angulo1-angulo2))<=angulo_max
    condicion_longitud=longitud1<=longitud_max | longitud2<=longitud_max
    
    ### Si se cumplen ambas condiciones se eliminan los ejes y se unen los dos puntos contiguos a v
    if (condicion_angulo & condicion_longitud){
      vertices_ejes_aux=vertices_ejes_aux[-ejes_eliminar,]
      vertices_ejes_aux=rbind(vertices_ejes_aux,c(puntos_unir[1],puntos_unir[2]))
    }
  }
  return(vertices_ejes_aux)
}

DibujarRed <- function(grafo_original,vertices_ejes_aux){
  indice=0
  lineas=list()
  for (i in c(1:nrow(vertices_ejes_aux))){
    indice=indice+1
    linea=Lines(list(Line(rbind(c(grafo_original$vertices$x[vertices_ejes_aux[i,1]],
                                  grafo_original$vertices$y[vertices_ejes_aux[i,1]]),
                                c(grafo_original$vertices$x[vertices_ejes_aux[i,2]],
                                  grafo_original$vertices$y[vertices_ejes_aux[i,2]])))),
                toString(indice))
    lineas[[indice]]=linea
  }
  lineas <- SpatialLines(lineas)
  return(lineas)
}

### GLOBAL VARIABLES (SOME ARE FINALLY NOT USED)

leyenda=0
valor_ok_global=0
contador_direccion=0
contador_direccion_dir=0
contador_direccion_point=0
spatial_lines_lin=NULL
spatial_lines_col=NULL
spatial_lines_global=NULL
spatial_lines_lin_dir=NULL
spatial_lines_col_dir=NULL
spatial_lines_global_dir=NULL
point_pattern_edited=NULL
ejes_eliminar=c()
vertices_unir=c()
puntos_unir=c()
puntos_unir_dos_nuevos=c()
Vertices=NULL
Vertices_dir=NULL
Vert=NULL
Vert_dir=NULL
Flow=c()
data_global=c()
Puntos_Actuales=c()
Puntos_Mover=c()
Marcas_Actuales=c()
Puntos_Marcas_Actuales=c()
object_global=NULL
state_popup_eventos_point=NULL
zoom_mapa=c()
center_mapa=c()
zoom_mapa_dir=c()
center_mapa_dir=c()
zoom_mapa_point=c()
center_mapa_point=c()
proj_edit_global=c()
proj_dir_global=c()

##################################################################################################################################
##########################################          SERVER          ##############################################################
##################################################################################################################################

shinyServer(function(input, output,session) {
  
  ### FUNCIONES CLICK
  
  data_of_click <- reactiveValues(clickedMarker=NULL)
  data_of_click_last <- reactiveValues(clickedMarkerLast=NULL)
  data_of_click_secondlast <- reactiveValues(clickedMarkerSecondLast=NULL)
  data_of_click_origin <- reactiveValues(clickedMarkerOrigin=NULL)
  data_of_click_origin_dir <- reactiveValues(clickedMarkerOrigin_dir=NULL)
  data_of_click_end <- reactiveValues(clickedMarkerEnd=NULL)
  data_of_click_end_dir <- reactiveValues(clickedMarkerEnd_dir=NULL)
  data_of_click_edge <- reactiveValues(clickedEdge=NULL)
  data_of_click_map <- reactiveValues(clickedPoint=NULL)
  data_of_click_map_secondlast <- reactiveValues(clickedPointSecondLast=NULL)
  data_of_click_map_dir <- reactiveValues(clickedPoint=NULL)
  data_of_click_remove_point <- reactiveValues(clickedMarkerRemove=NULL)
  data_of_click_new_point <- reactiveValues(clickedPointNew=NULL)
  
  data <- reactiveValues(clickedMarker=NULL)
  map = createLeafletMap(session, 'map')
  
##################################################################################################################################
########################################## NETWORDK EDITION ######################################################################
##################################################################################################################################
  
  
  output$mymap <- renderLeaflet({
    

      if (is.null(input$file)){
        return(NULL)
      } else{
        nombre=input$file
        print(nombre)
        print(nombre$datapath)
        spatial_lines<-readRDS(file=nombre$datapath)
        print(proj4string(spatial_lines))
        proj_edit_global<<-proj4string(spatial_lines)
      
        if (as.numeric(input$rebuildnet)==0){
          
        #### INPUT SpatialLines MODIFICATION
      
        spatial_lines_lin<<-as.linnet.SpatialLines(spatial_lines)
        vertices=vertices(spatial_lines_lin)
        spatial_lines <- spTransform(spatial_lines, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

        spatial_lines_lin<<-as.linnet.SpatialLines(spatial_lines)
        Vertices<<-vertices(spatial_lines_lin)
        Vert<<-cbind(spatial_lines_lin$from,spatial_lines_lin$to)
        spatial_lines_global<<-spatial_lines
        
        }
        
        ### GLOBAL SETTINGS MAP (DEFAULT)
        
        center_mapa<<-gCentroid(spatial_lines_global)
        zoom_mapa<<-13
        
        ### REBUILD IF BUTTON HAS BEEN CLICKED
        
        print(as.numeric(input$rebuildnet))
        if (as.numeric(input$rebuildnet)!=0){
          ReBuild_Net()
          ### CREAR NUEVO SpatialLines INCORPORANDO LOS EJES A ANYADIR
          lineas=list()
          indice=0
          for (i in c(1:spatial_lines_lin$lines$n)){
            indice=indice+1
            linea=Lines(list(Line(rbind(c(spatial_lines_lin$lines$ends$x0[i],
                                          spatial_lines_lin$lines$ends$y0[i]),
                                        c(spatial_lines_lin$lines$ends$x1[i],
                                          spatial_lines_lin$lines$ends$y1[i])))),
                        toString(indice))
            lineas[[indice]]=linea
          }
          lineas <- SpatialLines(lineas)
          
          spatial_lines_global<<-lineas

          Vertices<<-vertices(spatial_lines_lin)
          Vert<<-cbind(spatial_lines_lin$from,spatial_lines_lin$to)

        }
        
        ### SIMPLIFY IF BUTTON HAS BEEN CLICKED
        
        print(as.numeric(input$simplifynet))
        if (as.numeric(input$simplifynet)!=0){
          Simplify_Net()
          ### CREATE NEW SpatialLines
          lineas=list()
          indice=0
          print(class(spatial_lines_lin))
          for (i in c(1:spatial_lines_lin$lines$n)){
            indice=indice+1
            linea=Lines(list(Line(rbind(c(spatial_lines_lin$lines$ends$x0[i],
                                          spatial_lines_lin$lines$ends$y0[i]),
                                        c(spatial_lines_lin$lines$ends$x1[i],
                                          spatial_lines_lin$lines$ends$y1[i])))),
                        toString(indice))
            lineas[[indice]]=linea
          }
          lineas <- SpatialLines(lineas)
          # spatial_lines_global<<-SpatialLinesDataFrame(lineas, spatial_lines_global@data, match.ID = F)
          spatial_lines_global<<-lineas

          Vertices<<-vertices(spatial_lines_lin)
          Vert<<-cbind(spatial_lines_lin$from,spatial_lines_lin$to)

        }

        leaflet(spatial_lines_global) %>%
            addProviderTiles("Esri.WorldStreetMap") %>%
            setView(center_mapa$x, center_mapa$y, zoom = zoom_mapa) %>%
            addPolylines(color="black",
                         weight=input$edge_thickness)%>%
            addCircleMarkers(lng=Vertices$x,lat=Vertices$y,radius=input$vertex_radius,fillOpacity = 1,color="black")
      }
  })
  
  ##################################################################################################################################
  ########################################## NETWORK DIRECTION #####################################################################
  ##################################################################################################################################
  

  output$mymapdirection <- renderLeaflet({
    
    
    if (is.null(input$fileNetworkDirection)){
      return(NULL)
    } else{
      nombre=input$fileNetworkDirection
      print(nombre)
      print(nombre$datapath)
      spatial_lines_dir<-readRDS(file=nombre$datapath)
      print(proj4string(spatial_lines_dir))
      proj_dir_global<<-proj4string(spatial_lines_dir)
      
      #### INPUT SpatialLines MODIFICATION
      
      spatial_lines_lin_dir<<-as.linnet.SpatialLines(spatial_lines_dir)
      vertices_dir=vertices(spatial_lines_lin_dir)
      spatial_lines_dir <- spTransform(spatial_lines_dir, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
      spatial_lines_lin_dir<<-as.linnet.SpatialLines(spatial_lines_dir)
      Vertices_dir<<-vertices(spatial_lines_lin_dir)
      Vert_dir<<-cbind(spatial_lines_lin_dir$from,spatial_lines_lin_dir$to)
        
      df=data.frame(Vert_dir,rep(0,spatial_lines_lin_dir$lines$n))
      df=data.frame(Vert_dir,rep(0,length(spatial_lines_dir)))
      colnames(df)=c("V1","V2","Dir")

      Flow<<-rep(0,spatial_lines_lin_dir$lines$n)

      ### CREATE SpatialLinesDataFrame IF INPUT HAS NO DATA.FRAME ATTACHED
      ### IF DATA.FRAME ATTACHED BUT NO Dir OR V1 OR V2, THEN ADD

      if (class(spatial_lines_dir)=="SpatialLines"){
        spatial_lines_dir <- SpatialLinesDataFrame(spatial_lines_dir, df, match.ID = F)
        data_global_dir<<-df
        print(class(spatial_lines_dir))
      } else{
        if (!"V1" %in% names(spatial_lines_dir@data) | !"V2" %in% names(spatial_lines_dir@data) |
            !"Dir" %in% names(spatial_lines_dir@data)){
          spatial_lines_dir <- SpatialLinesDataFrame(spatial_lines_dir, data.frame(spatial_lines_dir@data,df), 
                                                     match.ID = F)
          data_global_dir<<-data.frame(spatial_lines_dir@data,df)
        }
      }
        
      ### EDGE INFO WHEN CLICK
        
      spatial_lines_global_dir<<-spatial_lines_dir
      
      directions_edges=as.numeric(spatial_lines_global_dir@data$Dir)
      directed_edges_1=which((directions_edges==1)==T)
      directed_edges_menos1=which((directions_edges==-1)==T)
      directed_edges_2=which((directions_edges==2)==T)
      directed_edges_1=c(directed_edges_1,directed_edges_2)
      directed_edges_menos1=c(directed_edges_menos1,directed_edges_2)
   
      ### GLOBAL SETTINGS 
      
      center_mapa_dir=gCentroid(spatial_lines_global_dir)
   
      if (length(directed_edges_1)==0 & length(directed_edges_menos1)==0 & length(directed_edges_2)==0){
        
        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")
      } else if (length(directed_edges_1)>0 & length(directed_edges_menos1)==0 & length(directed_edges_2)==0){
        
        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_1,2]], 
                   color = "#0078ff", maxThickness = 1.6, dir=1)
      } else if (length(directed_edges_1)==0 & length(directed_edges_menos1)>0 & length(directed_edges_2)==0){

        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_menos1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_menos1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_menos1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_menos1,2]], 
                   color = "#0078ff", maxThickness = 1.6, dir=-1)
      } else if (length(directed_edges_1)==0 & length(directed_edges_menos1)==0 & length(directed_edges_2)>0){
        
        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]], 
                   color = "#e1009a", maxThickness = 1.6, dir=0) %>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]], 
                   color = "#e1009a", maxThickness = 1.6, dir=0)
      } else if (length(directed_edges_1)>0 & length(directed_edges_menos1)>0 & length(directed_edges_2)==0){

        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_1,2]],
                   color = "#0078ff", maxThickness = 1.6, dir=1)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_menos1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_menos1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_menos1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_menos1,2]],
                   color = "#0078ff", maxThickness = 1.6, dir=-1)
      } else if (length(directed_edges_1)==0 & length(directed_edges_menos1)>0 & length(directed_edges_2)>0){
        
        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_menos1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_menos1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_menos1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_menos1,2]],
                   color = "#0078ff", maxThickness = 1.6, dir=-1)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]],
                   color = "#e1009a", maxThickness = 1.6, dir=0)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]],
                   color = "#e1009a", maxThickness = 1.6, dir=0)
      } else if (length(directed_edges_1)>0 & length(directed_edges_menos1)==0 & length(directed_edges_2)>0){
        
        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_1,2]],
                   color = "#0078ff", maxThickness = 1.6, dir=-1)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]],
                   color = "#e1009a", maxThickness = 1.6, dir=0)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]],
                   color = "#e1009a", maxThickness = 1.6, dir=0)
      } else if (length(directed_edges_1)>0 & length(directed_edges_menos1)>0 & length(directed_edges_2)>0){
        leaflet(spatial_lines_global_dir) %>%
          addTiles() %>%
          setView(center_mapa_dir$x, center_mapa_dir$y, zoom = 13) %>%
          addPolylines(color="#545454",
                       weight=5)%>%
          addCircleMarkers(lng=Vertices_dir$x,lat=Vertices_dir$y,radius=8,fillOpacity = 1,color="black")%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_1,2]],
                   color = "#0078ff", maxThickness = 1.6, dir=1)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_menos1,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_menos1,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_menos1,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_menos1,2]],
                   color = "#0078ff", maxThickness = 1.6, dir=-1)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]],
                   color = "#e1009a", maxThickness = 1.6, dir=0)%>%
          addFlows(lng0=Vertices_dir$x[Vert_dir[directed_edges_2,1]],lat0=Vertices_dir$y[Vert_dir[directed_edges_2,1]],
                   lng1=Vertices_dir$x[Vert_dir[directed_edges_2,2]],lat1=Vertices_dir$y[Vert_dir[directed_edges_2,2]],
                   color = "#e1009a", maxThickness = 1.6, dir=0)
      }
    }
  })
  
  ##################################################################################################################################
  ########################################## AUXILIARY FOR EDITION AND DIRECTION ###################################################
  ##################################################################################################################################
  

  observeEvent(input$mymap_shape_click, {
    data_of_click_edge$clickedEdge <- input$mymap_shape_click
  })
  
  ### CONTROL MARKER ORIGIN AND END
  observeEvent(input$mymap_marker_click, {
    data_of_click_secondlast$clickedMarkerSecondLast <- data_of_click_last$clickedMarkerLast
    data_of_click_last$clickedMarkerLast <- input$mymap_marker_click
    if (contador_direccion%%2==0){
      click <- input$mymap_marker_click
      data_of_click_origin$clickedMarkerOrigin <- input$mymap_marker_click
      print(data_of_click_origin$clickedMarkerOrigin$lng)
      print(data_of_click_origin$clickedMarkerOrigin$lat)
      contador_direccion<<-contador_direccion+1
      # print(paste0("contadorORIGIN",contador_direccion))
    } else {
      click <- input$mymap_marker_click
      data_of_click_end$clickedMarkerEnd <- input$mymap_marker_click
      print(data_of_click_end$clickedMarkerEnd$lng)
      print(data_of_click_end$clickedMarkerEnd$lat)
      contador_direccion<<-contador_direccion+1
      # print(paste0("contadorEND",contador_direccion))
    }
  })
  
  ### CONTROL MARKER ORIGIN AND END DIRECTION
  observeEvent(input$mymapdirection_marker_click, {
    if (contador_direccion_dir%%2==0){
      click <- input$mymapdirection_marker_click
      data_of_click_origin_dir$clickedMarkerOrigin_dir <- input$mymapdirection_marker_click
      contador_direccion_dir<<-contador_direccion_dir+1
    } else {
      click <- input$mymapdirection_marker_click
      data_of_click_end_dir$clickedMarkerEnd_dir <- input$mymapdirection_marker_click
      contador_direccion_dir<<-contador_direccion_dir+1
    }
  })
  
  observeEvent(input$mymap_click,{
    ### ONLY IF IT IS NOT A MARKER
    Point=list(x=input$mymap_click$lng,y=input$mymap_click$lat)
    V_Point=BuscarVertice(Vertices,Point,10^(-6))
    print(paste0("V_Point",V_Point))
    if (length(V_Point)==0){
      if (!is.null(data_of_click_map$clickedPoint)){
        data_of_click_map_secondlast$clickedPointSecondLast <- data_of_click_map$clickedPoint
      }
      data_of_click_map$clickedPoint <- input$mymap_click
    }
    # print(data_of_click_map$clickedPoint$lng)
    # print(data_of_click_map$clickedPoint$lat)
    # print(data_of_click_map_secondlast$clickedPointSecondLast$lng)
    # print(data_of_click_map_secondlast$clickedPointSecondLast$lat)
  })
  
  ### CLEAN ALL STORED CLICKS WHEN ACTION IS CHANGED
  observeEvent(input$action,{
    data_of_click_origin$clickedMarkerOrigin <- NULL
    data_of_click_end$clickedMarkerEnd <- NULL
    data_of_click_last$clickedMarkerLast <- NULL
    data_of_click_edge$clickedEdge <- NULL
    data_of_click_map$clickedPoint <- NULL
    data_of_click_map_secondlast$clickedPointSecondLast <- NULL
    # print(data_of_click_origin$clickedMarkerOrigin$lng)
    # print(data_of_click_origin$clickedMarkerOrigin$lat)
    # print(data_of_click_end$clickedMarkerEnd$lng)
    # print(data_of_click_end$clickedMarkerEnd$lat)
  })
  
  observeEvent(input$action_flow,{
    data_of_click_origin_dir$clickedMarkerOrigin_dir <- NULL
    data_of_click_end_dir$clickedMarkerEnd_dir <- NULL
    data_of_click_map_dir$clickedPoint <- NULL
    contador_direccion_dir <<- 0
    # buscar_eje=NA
  })
  
  ### ADD ARROWS TO THE NETWORK
  observe({
 
    if (!is.null(data_of_click_origin_dir$clickedMarkerOrigin_dir) & !is.null(data_of_click_end_dir$clickedMarkerEnd_dir) & contador_direccion_dir%%2==0){
      
      fromPoint=list(x=data_of_click_origin_dir$clickedMarkerOrigin_dir$lng,y=data_of_click_origin_dir$clickedMarkerOrigin_dir$lat)
      toPoint=list(x=data_of_click_end_dir$clickedMarkerEnd_dir$lng,y=data_of_click_end_dir$clickedMarkerEnd_dir$lat)
      
      ### FIND ORIGIN AND END
  
      if (input$action_flow=="add.flow"){
        
        ### FIND VERTEX AND MODIFY DATA.FRAME
        
        V_origen=BuscarVertice(Vertices_dir,fromPoint,10^(-6))
        V_fin=BuscarVertice(Vertices_dir,toPoint,10^(-6))
        # print(V_origen)
        # print(V_fin)

        buscar_eje=row.match(c(V_origen,V_fin),Vert_dir)
        # print(buscar_eje)

        if (!is.na(buscar_eje)){
          if (spatial_lines_global_dir@data$Dir[buscar_eje]==0){
            spatial_lines_global_dir@data$Dir[buscar_eje] <<- 1
          } 
          else if (spatial_lines_global_dir@data$Dir[buscar_eje]==-1){
            spatial_lines_global_dir@data$Dir[buscar_eje] <<- 2
          }
          
          # print(spatial_lines_global_dir@data$Dir[buscar_eje])
        } else{
          buscar_eje=row.match(c(V_fin,V_origen),Vert_dir)
          if (!is.na(buscar_eje)){
            # print(buscar_eje)
            # print(spatial_lines_lin)
            if (spatial_lines_global_dir@data$Dir[buscar_eje]==0){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- -1
            }
            else if (spatial_lines_global_dir@data$Dir[buscar_eje]==1){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 2
            }
            # print(spatial_lines_global_dir@data$Dir[buscar_eje])
          }
        }
        
        leafletProxy("mymapdirection")%>%addFlows(lng0=fromPoint$x,lat0=fromPoint$y,lng1=toPoint$x,lat1=toPoint$y, color = "#0078ff",
                                           maxThickness = 1.6, dir=1)
      } 
      if (input$action_flow=="remove.flow"){
        
        ### FIND VERTEX AND MODIFY DATA.FRAME
        
        V_origen=BuscarVertice(Vertices_dir,fromPoint,10^(-6))
        V_fin=BuscarVertice(Vertices_dir,toPoint,10^(-6))
        # print(paste0("V_origen",V_origen))
        # print(paste0("V_fin",V_fin))

        buscar_eje=row.match(c(V_origen,V_fin),Vert_dir)
        # print(paste0("buscar_eje",buscar_eje))
        
        if (!is.na(buscar_eje)){
          if (spatial_lines_global_dir@data$Dir[buscar_eje]==1){
            spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
            # print(spatial_lines_global_dir@data$Dir[buscar_eje])
          } 
          else if (spatial_lines_global_dir@data$Dir[buscar_eje]==-1){
            spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
            # print(spatial_lines_global_dir@data$Dir[buscar_eje])
          } 
          else if (spatial_lines_global_dir@data$Dir[buscar_eje]==2){
            spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
            # print(spatial_lines_global_dir@data$Dir[buscar_eje])
          }
        } else{
          buscar_eje=row.match(c(V_fin,V_origen),Vert_dir)
          aux_V_origen=V_origen
          V_origen=V_fin
          V_fin=aux_V_origen
          aux_fromPoint_x=fromPoint$x
          aux_fromPoint_y=fromPoint$y
          fromPoint$x=toPoint$x
          fromPoint$y=toPoint$y
          toPoint$x=aux_fromPoint_x
          toPoint$y=aux_fromPoint_y
          # print(paste0("V_origen",V_origen))
          # print(paste0("V_fin",V_fin))
          # print(paste0("buscar_eje",buscar_eje))
          # print(spatial_lines_lin_dir)
          if (!is.na(buscar_eje)){
            if (spatial_lines_global_dir@data$Dir[buscar_eje]==-1){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
              # print(spatial_lines_global_dir@data$Dir[buscar_eje])
            }
            else if (spatial_lines_global_dir@data$Dir[buscar_eje]==1){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
              # print(spatial_lines_global_dir@data$Dir[buscar_eje])
            }
            else if (spatial_lines_global_dir@data$Dir[buscar_eje]==2){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
              # print(spatial_lines_global_dir@data$Dir[buscar_eje])
            }
          }
        }
        
        leafletProxy("mymapdirection")%>%
          addFlows(lng0=fromPoint$x,lat0=fromPoint$y,lng1=toPoint$x,lat1=toPoint$y,
                   color = "#545454",maxThickness = 1.7, dir=0)
        
        # leafletProxy("mymapdirection")%>%addFlows(lng0=Vertices_dir$x[V_origen],lat0=Vertices_dir$y[V_origen],
        #                                           lng1=Vertices_dir$x[V_fin],lat1=Vertices_dir$y[V_fin],
        #                                           color = "black",
        #                                           maxThickness = 0, dir=1)
      

        # leafletProxy("mymapdirection")%>%addFlows(lng0=Vertices_dir$x[Vert_dir[,]],
        #                                           lat0=Vertices_dir$y[Vert_dir[,]],
        #                                           lng1=Vertices_dir$x[Vert_dir[,]],
        #                                           lat1=Vertices_dir$y[Vert_dir[,]],
        #                                           color = "black",
        #                                           maxThickness = 0, dir=0)
        
        
      }
      
      if (input$action_flow=="add.longflow"){
        V_origen=BuscarVertice(Vertices_dir,fromPoint,10^(-6))
        V_fin=BuscarVertice(Vertices_dir,toPoint,10^(-6))
        print(V_origen)
        print(V_fin)
        flow_vertex=shortestpath(spatial_lines_lin_dir,V_origen,V_fin)
        print(flow_vertex)
        for (j in c(1:(length(flow_vertex)-1))){
          
          V_origen=flow_vertex[j]
          V_fin=flow_vertex[j+1]
          
          buscar_eje=row.match(c(V_origen,V_fin),Vert_dir)
          print(buscar_eje)
          
          if (!is.na(buscar_eje)){
            if (spatial_lines_global_dir@data$Dir[buscar_eje]==0){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 1
            } 
            else if (spatial_lines_global_dir@data$Dir[buscar_eje]==-1){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 2
            }
            print(spatial_lines_global_dir@data$Dir[buscar_eje])
          } else{
            buscar_eje=row.match(c(V_fin,V_origen),Vert_dir)
            print(buscar_eje)
            print(spatial_lines_lin_dir)
            if (!is.na(buscar_eje)){
              if (spatial_lines_global_dir@data$Dir[buscar_eje]==0){
                spatial_lines_global_dir@data$Dir[buscar_eje] <<- -1
              }
              else if (spatial_lines_global_dir@data$Dir[buscar_eje]==1){
                spatial_lines_global_dir@data$Dir[buscar_eje] <<- 2
              }
              print(spatial_lines_global_dir@data$Dir[buscar_eje])
            }
          }
        }
        ### ARROW ADDITION
        for (j in c(1:(length(flow_vertex)-1))){
          leafletProxy("mymapdirection")%>%addFlows(lng0=Vertices_dir$x[flow_vertex[j]],lat0=Vertices_dir$y[flow_vertex[j]],
                                           lng1=Vertices_dir$x[flow_vertex[j+1]],lat1=Vertices_dir$y[flow_vertex[j+1]], color = "#0078ff",
                                           maxThickness = 1.6, dir=1)
        }
      }
      
      if (input$action_flow=="remove.longflow"){
        V_origen=BuscarVertice(Vertices_dir,fromPoint,10^(-6))
        V_fin=BuscarVertice(Vertices_dir,toPoint,10^(-6))
        print(V_origen)
        print(V_fin)
        flow_vertex=shortestpath(spatial_lines_lin_dir,V_origen,V_fin)
        print(flow_vertex)
        for (j in c(1:(length(flow_vertex)-1))){
          
          V_origen=flow_vertex[j]
          V_fin=flow_vertex[j+1]
          
          buscar_eje=row.match(c(V_origen,V_fin),Vert_dir)
          print(buscar_eje)
          
          if (!is.na(buscar_eje)){
            if (spatial_lines_global_dir@data$Dir[buscar_eje]==1){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
            } 
            else if (spatial_lines_global_dir@data$Dir[buscar_eje]==2){
              spatial_lines_global_dir@data$Dir[buscar_eje] <<- -1
            }
            print(spatial_lines_global_dir@data$Dir[buscar_eje])
          } else{
            buscar_eje=row.match(c(V_fin,V_origen),Vert_dir)
            print(buscar_eje)
            print(spatial_lines_lin_dir)
            if (!is.na(buscar_eje)){
              if (spatial_lines_global_dir@data$Dir[buscar_eje]==-1){
                spatial_lines_global_dir@data$Dir[buscar_eje] <<- 0
              }
              else if (spatial_lines_global_dir@data$Dir[buscar_eje]==2){
                spatial_lines_global_dir@data$Dir[buscar_eje] <<- 1
              }
              print(spatial_lines_global_dir@data$Dir[buscar_eje])
            }
          }
        }
        ### ARROW REMOVAL
        for (j in c(1:(length(flow_vertex)-1))){
          leafletProxy("mymapdirection")%>%addFlows(lng0=Vertices_dir$x[flow_vertex[j]],lat0=Vertices_dir$y[flow_vertex[j]],
                                           lng1=Vertices_dir$x[flow_vertex[j+1]],lat1=Vertices_dir$y[flow_vertex[j+1]], color = "#545454",
                                           maxThickness = 1.7, dir=0)
        }
      }
    }
  })
  
  ### REMOVE EDGE NETWORK
  observe({
    if (input$action=="remove.edge" & !is.null(data_of_click_edge$clickedEdge)){
      eje_eliminar=as.numeric(DetectarEjeGrafo(spatial_lines_lin,list(x=data_of_click_edge$clickedEdge$lng,
                                                           y=data_of_click_edge$clickedEdge$lat)))
      x_eliminar=c(Vertices$x[Vert[eje_eliminar,1]],Vertices$x[Vert[eje_eliminar,2]])
      y_eliminar=c(Vertices$y[Vert[eje_eliminar,1]],Vertices$y[Vert[eje_eliminar,2]])
      leafletProxy("mymap")%>%addPolylines(lng=x_eliminar,lat=y_eliminar, color = "red")
      ejes_eliminar<<-c(ejes_eliminar,eje_eliminar)
    }
  })
  
  ### ADD EDGE NETWORK (JOIN VERTEX)
  observe({
    if (input$action=="add.edge" & !is.null(data_of_click_origin$clickedMarkerOrigin) & !is.null(data_of_click_end$clickedMarkerEnd)){
      x_unir=c(data_of_click_origin$clickedMarkerOrigin$lng,data_of_click_end$clickedMarkerEnd$lng)
      y_unir=c(data_of_click_origin$clickedMarkerOrigin$lat,data_of_click_end$clickedMarkerEnd$lat)
      V_origen=BuscarVertice(Vertices,list(x=x_unir[1],y=y_unir[1]),10^(-6))
      V_fin=BuscarVertice(Vertices,list(x=x_unir[2],y=y_unir[2]),10^(-6))
      leafletProxy("mymap")%>%addPolylines(lng=x_unir,lat=y_unir, color = "green")
      vertices_unir<<-rbind(vertices_unir,c(V_origen,V_fin))
      data_of_click_origin$clickedMarkerOrigin <- NULL
      data_of_click_end$clickedMarkerEnd <- NULL
    }
  })
  
  ### ADD POINT + EDGE
  observe({
    if (input$action=="add.point" & !is.null(data_of_click_last$clickedMarkerLast) & !is.null(data_of_click_map$clickedPoint)){
      x_unir=c(data_of_click_last$clickedMarkerLast$lng,data_of_click_map$clickedPoint$lng)
      y_unir=c(data_of_click_last$clickedMarkerLast$lat,data_of_click_map$clickedPoint$lat)
      ### FIND VERTEX
      V_join=BuscarVertice(Vertices,list(x=x_unir[1],y=y_unir[1]),10^(-6))
      leafletProxy("mymap")%>%addPolylines(lng=x_unir,lat=y_unir, color = "green")
      puntos_unir<<-rbind(puntos_unir,c(V_join,x_unir[2],y_unir[2]))
      data_of_click_map$clickedPoint=NULL
      data_of_click_last$clickedMarkerLast=NULL
    }
  })
  
  ### ADD 2 POINTS + EDGE
  observe({
    if (input$action=="add.twopoints" & !is.null(data_of_click_map_secondlast$clickedPointSecondLast) & !is.null(data_of_click_map$clickedPoint)){
      x_unir=c(data_of_click_map_secondlast$clickedPointSecondLast$lng,data_of_click_map$clickedPoint$lng)
      y_unir=c(data_of_click_map_secondlast$clickedPointSecondLast$lat,data_of_click_map$clickedPoint$lat)
      leafletProxy("mymap")%>%addPolylines(lng=x_unir,lat=y_unir, color = "green")
      puntos_unir_dos_nuevos<<-rbind(puntos_unir_dos_nuevos,c(x_unir[1],y_unir[1],x_unir[2],y_unir[2]))
      data_of_click_map$clickedPoint=NULL
      data_of_click_map_secondlast$clickedPointSecondLast=NULL
    }
  })
  
  ### NETWORK SIMPLIFICATION
  Simplify_Net <- eventReactive(input$simplifynet, {
    if (!is.na(as.numeric(input$angle_max)) & !is.na(as.numeric(input$length_max))){
      aux_linnet=spatial_lines_lin
      aux_simplify=SimplificarRed(aux_linnet,input$angle_max,input$length_max,Vert)
      aux_simplify=DibujarRed(aux_linnet,aux_simplify)
      aux_linnet=as.linnet(aux_simplify)
      spatial_lines_lin<<-aux_linnet
    }
    return(spatial_lines_lin)
  })
  
  ### REBUILDING THE NETWORK
  
  ReBuild_Net <- eventReactive(input$rebuildnet, {
    
    ### UPDATE ZOOM LEVEL
    zoom_mapa<<-input$mymap_zoom
    bounds_mapa=input$mymap_bounds
    center_mapa<<-list(x=mean(c(bounds_mapa$east, bounds_mapa$west)),y=mean(c(bounds_mapa$north, bounds_mapa$south)))
    aux_linnet=spatial_lines_lin
  
    ### FIRST, JOIN VERTEX
    
    if (!is.null(vertices_unir)){
      aux_sp=list()
      n_ejes=aux_linnet$lines$n
      for (i in c(1:(n_ejes+nrow(vertices_unir)))){
        if (i<=n_ejes){
          linea=Lines(list(Line(rbind(c(aux_linnet$lines$ends$x0[i],aux_linnet$lines$ends$y0[i]),
                                      c(aux_linnet$lines$ends$x1[i],aux_linnet$lines$ends$y1[i])))), toString(i))
          aux_sp[[i]]=linea
        } else {
          j=i-n_ejes
          linea=Lines(list(Line(rbind(c(aux_linnet$vertices$x[vertices_unir[j,1]],aux_linnet$vertices$y[vertices_unir[j,1]]),
                                      c(aux_linnet$vertices$x[vertices_unir[j,2]],aux_linnet$vertices$y[vertices_unir[j,2]])))), toString(i))
          aux_sp[[i]]=linea
        }
        
      }
      aux_sp <- SpatialLines(aux_sp)
      aux_linnet=as.linnet(aux_sp)
      ### CLEAN VECTOR
      vertices_unir<<-c()
    }
    
    ### SECOND, ADD A NEW POINT
    
    if (!is.null(puntos_unir)){
      aux_sp=list()
      n_ejes=aux_linnet$lines$n
      for (i in c(1:(n_ejes+nrow(puntos_unir)))){
        if (i<=n_ejes){
          linea=Lines(list(Line(rbind(c(aux_linnet$lines$ends$x0[i],aux_linnet$lines$ends$y0[i]),
                                      c(aux_linnet$lines$ends$x1[i],aux_linnet$lines$ends$y1[i])))), toString(i))
          aux_sp[[i]]=linea
        } else {
          j=i-n_ejes
          linea=Lines(list(Line(rbind(c(aux_linnet$vertices$x[puntos_unir[j,1]],aux_linnet$vertices$y[puntos_unir[j,1]]),
                                      c(puntos_unir[j,2],puntos_unir[j,3])))), toString(i))
          aux_sp[[i]]=linea
        }
        
      }
      aux_sp <- SpatialLines(aux_sp)
      aux_linnet=as.linnet(aux_sp)
      ### CLEAN VECTOR
      puntos_unir<<-c()
    }
    
    ### THIRD, ADD TWO NEW POINTS
    
    if (!is.null(puntos_unir_dos_nuevos)){
      aux_sp=list()
      n_ejes=aux_linnet$lines$n
      for (i in c(1:(n_ejes+nrow(puntos_unir_dos_nuevos)))){
        if (i<=n_ejes){
          linea=Lines(list(Line(rbind(c(aux_linnet$lines$ends$x0[i],aux_linnet$lines$ends$y0[i]),
                                      c(aux_linnet$lines$ends$x1[i],aux_linnet$lines$ends$y1[i])))), toString(i))
          aux_sp[[i]]=linea
        } else {
          j=i-n_ejes
          linea=Lines(list(Line(rbind(c(puntos_unir_dos_nuevos[j,1],puntos_unir_dos_nuevos[j,2]),
                                      c(puntos_unir_dos_nuevos[j,3],puntos_unir_dos_nuevos[j,4])))), toString(i))
          aux_sp[[i]]=linea
        }
        
      }
      aux_sp <- SpatialLines(aux_sp)
      aux_linnet=as.linnet(aux_sp)
      ### CLEAN VECTOR
      puntos_unir_dos_nuevos<<-c()
    }
    
    ### LAST, REMOVE EDGES
    
    if (!is.null(ejes_eliminar)){
      edges_retain=c(1:aux_linnet$lines$n)
      print(edges_retain)
      print(ejes_eliminar)
      edges_retain=edges_retain[-ejes_eliminar]
      aux_linnet=thinNetwork(aux_linnet,retainedges = edges_retain)
      ### CLEAN VECTOR
      ejes_eliminar<<-c()
      Vertices<<-vertices(spatial_lines_lin)
      Vert<<-cbind(spatial_lines_lin$from,spatial_lines_lin$to)
    }
    ### RETURN
    spatial_lines_lin<<-aux_linnet
    return(spatial_lines_lin)
  })
  
  output$downloadnet <- downloadHandler(
    filename = function() {
      paste("EditedRoadNetwork ",gsub(":","_",date()),".rds", sep = "")
    },
    content = function(file) {
      proj4string(spatial_lines_global)=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      spatial_lines_download <- spTransform(spatial_lines_global, CRS(proj_edit_global))
      saveRDS(spatial_lines_download, file)
    }
  )
  
  output$downloadnetNetworkDirection <- downloadHandler(
    filename = function() {
      paste("DirectedRoadNetwork ",gsub(":","_",date()),".rds", sep = "")
    },
    content = function(file) {
      proj4string(spatial_lines_global_dir)=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      spatial_lines_download_dir <- spTransform(spatial_lines_global_dir, CRS(proj_dir_global))
      saveRDS(spatial_lines_download_dir, file)
    }
  )
 
  ##################################################################################################################################
  ########################################## POINT PATTERN EDITION/ VISUALIZATION #########################################################################
  ##################################################################################################################################
  
  output$mymapPointPatternEdition <- renderLeaflet({
    
    if (is.null(input$filePointPattern_lppx)){
      return(NULL)
    } else{
      nombre=input$filePointPattern_lppx
      print(nombre)
      print(nombre$datapath)

      object<-readRDS(file=nombre$datapath)
      object_global<<-object
      
      print(class(object))

      lineas=list()
      indice=0
      for (i in c(1:object$domain$lines$n)){
        indice=indice+1
        linea=Lines(list(Line(rbind(c(object$domain$lines$ends$x0[i],
                                      object$domain$lines$ends$y0[i]),
                                    c(object$domain$lines$ends$x1[i],
                                      object$domain$lines$ends$y1[i])))),
                    toString(indice))
        lineas[[indice]]=linea
      }
      lineas <- SpatialLines(lineas)
      spatial_lines_global_point_pattern<<-lineas

      #### PROJECTION FIXING AND CONVERSION AND POINT PATTERN EXTRACTION
      
      proj4string(spatial_lines_global_point_pattern)=CRS(paste0("+proj=utm +zone=",input$utm_select," ellps=WGS84"))
      spatial_lines_global_point_pattern <- spTransform(spatial_lines_global_point_pattern, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      object_domain=as.linnet(spatial_lines_global_point_pattern)
      Vertices=vertices(object_domain)
        
      Puntos=cbind(object$data$x,object$data$y)
      Puntos=UTM2LONLAT(Puntos,zone=input$utm_select)

      Puntos_Actuales<<-cbind(Puntos[,1],Puntos[,2])
      print(paste0("PUNTOS ACTUALES ",Puntos_Actuales[1108,]))
      
      if (!is.null(marks(object))){
        Marcas_Actuales<<-data.frame(marks(object))
        Puntos_Marcas_Actuales<<-cbind(Puntos_Actuales,Marcas_Actuales)
      }
        
      if (!is.null(marks(object))){
        Puntos=data.frame(Puntos,marks(object))

        ### DETECT NUMERICAL MARKS
          
        numerical_marks=c()
        for (l in c(1:ncol(Puntos))){
          if (class(Puntos[,l])=="numeric" | class(Puntos[,l])=="integer"){
            numerical_marks=c(numerical_marks,l)
          }
        }
        for (l in c(1:length(numerical_marks))){
          if (numerical_marks[l]>=3){
            Puntos[,numerical_marks[l]]=round(Puntos[,numerical_marks[l]],digits=2)
          }
        }

        ### POPUP WHEN EVENT CLICKED
          
        if (ncol(Puntos_Marcas_Actuales)==3){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3])
        }
        if (ncol(Puntos_Marcas_Actuales)==4){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                       '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,3])
        }
        if (ncol(Puntos_Marcas_Actuales)==5){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5])
        }
        if (ncol(Puntos_Marcas_Actuales)==6){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[6],": ",'</strong>',Puntos_Marcas_Actuales[,6])
        }
        if (ncol(Puntos_Marcas_Actuales)==7){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[6],": ",'</strong>',Puntos_Marcas_Actuales[,6],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[7],": ",'</strong>',Puntos_Marcas_Actuales[,7])
        }
        if (ncol(Puntos_Marcas_Actuales)==8){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[6],": ",'</strong>',Puntos_Marcas_Actuales[,6],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[7],": ",'</strong>',Puntos_Marcas_Actuales[,7],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[8],": ",'</strong>',Puntos_Marcas_Actuales[,8])
        }
        if (ncol(Puntos_Marcas_Actuales)==9){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[6],": ",'</strong>',Puntos_Marcas_Actuales[,6],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[7],": ",'</strong>',Puntos_Marcas_Actuales[,7],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[8],": ",'</strong>',Puntos_Marcas_Actuales[,8],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[9],": ",'</strong>',Puntos_Marcas_Actuales[,9])
        }
        if (ncol(Puntos_Marcas_Actuales)==10){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[6],": ",'</strong>',Puntos_Marcas_Actuales[,6],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[7],": ",'</strong>',Puntos_Marcas_Actuales[,7],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[8],": ",'</strong>',Puntos_Marcas_Actuales[,8],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[9],": ",'</strong>',Puntos_Marcas_Actuales[,9],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[10],": ",'</strong>',Puntos_Marcas_Actuales[,10])
        }
        if (ncol(Puntos_Marcas_Actuales)>=11){
          state_popup_eventos_point <<- paste('<strong>',colnames(Puntos_Marcas_Actuales)[3],": ",'</strong>',Puntos_Marcas_Actuales[,3],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[4],": ",'</strong>',Puntos_Marcas_Actuales[,4],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[5],": ",'</strong>',Puntos_Marcas_Actuales[,5],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[6],": ",'</strong>',Puntos_Marcas_Actuales[,6],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[7],": ",'</strong>',Puntos_Marcas_Actuales[,7],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[8],": ",'</strong>',Puntos_Marcas_Actuales[,8],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[9],": ",'</strong>',Puntos_Marcas_Actuales[,9],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[10],": ",'</strong>',Puntos_Marcas_Actuales[,10],"<dd>",
                                         '<strong>',colnames(Puntos_Marcas_Actuales)[11],": ",'</strong>',Puntos_Marcas_Actuales[,11])
        }
      }
      
      ### GLOBAL SETTINGS MAP (DEFAULT)
      
      center_mapa_point<<-gCentroid(spatial_lines_global_point_pattern)
      zoom_mapa_point<<-13

      if (as.numeric(input$rebuildpoint)!=0){
        ### REBUILD POINT PATTERN
        ReBuild_PointPattern()
        
        object=point_pattern_edited
        object_global<<-object
        
        lineas=list()
        indice=0
        for (i in c(1:object$domain$lines$n)){
          indice=indice+1
          linea=Lines(list(Line(rbind(c(object$domain$lines$ends$x0[i],
                                        object$domain$lines$ends$y0[i]),
                                      c(object$domain$lines$ends$x1[i],
                                        object$domain$lines$ends$y1[i])))),
                      toString(indice))
          lineas[[indice]]=linea
        }
        lineas <- SpatialLines(lineas)
        spatial_lines_global_point_pattern<<-lineas
      }

      if (input$clusterevents_edition=="Yes"){
        leaflet(spatial_lines_global_point_pattern) %>%
          addProviderTiles("Esri.WorldStreetMap") %>%
          setView(center_mapa_point$x, center_mapa_point$y, zoom = zoom_mapa_point) %>%
          addPolylines(color="black",
                       weight=5)%>%
          addCircleMarkers(lng=Puntos_Actuales[,1],lat=Puntos_Actuales[,2],radius=4,fillOpacity = 1,color="red",
                           clusterOptions = markerClusterOptions(freezeAtZoom = F),
                           popup = state_popup_eventos_point,
                           popupOptions = popupOptions(maxWidth ="100%", closeOnClick = F))
      } else{
        leaflet(spatial_lines_global_point_pattern) %>%
          addProviderTiles("Esri.WorldStreetMap") %>% 
          setView(center_mapa_point$x, center_mapa_point$y, zoom = zoom_mapa_point) %>%
          addPolylines(color= "black",
                       weight=5)%>%
          addCircleMarkers(lng=Puntos_Actuales[,1],lat=Puntos_Actuales[,2],radius=4,fillOpacity = 1,color="red",
                           popup = state_popup_eventos_point,
                           popupOptions = popupOptions(maxWidth ="100%", closeOnClick = F))
      }
    }
  })
  
  output$downloadpointpattern <- downloadHandler(
    filename = function() {
      paste("EditedPointPattern ",gsub(":","_",date()),".rds", sep = "")
    },
    content = function(file) {
      incidentes_ppp=ppp(x=object_global$data$x,y=object_global$data$y, 
                           c(object_global$domain$window$xrange[1]-1000,object_global$domain$window$xrange[2]+1000),
                           c(object_global$domain$window$yrange[1]-1000,object_global$domain$window$yrange[2]+1000))
      point_pattern_download<<-lpp(incidentes_ppp,object_global$domain)
      marks(point_pattern_download)=data.frame(marks(object_global))
      saveRDS(point_pattern_download, file)
    }
  )
  
  ### CONTROL MARKER TO REMOVE POINT PATTERN
  observeEvent(input$mymapPointPatternEdition_marker_click, {
      click <- input$mymapPointPatternEdition_marker_click
      print(click)
      data_of_click_remove_point$clickedMarkerRemove <- input$mymapPointPatternEdition_marker_click
      print(paste("data_of_click_remove_point",data_of_click_remove_point$clickedMarkerRemove))
      contador_direccion_point<<-contador_direccion_point+1
  })
  
  ### CONTROL POINT TO ADD TO POINT PATTERN
  observeEvent(input$mymapPointPatternEdition_click, {
      print("ENTRA MAP CLICK")
      click <- input$mymapPointPatternEdition_click
      print(click)
      data_of_click_new_point$clickedPointNew <- input$mymapPointPatternEdition_click
      print(paste("data_of_click_new_point$clickedPointNew",data_of_click_new_point$clickedPointNew))
      contador_direccion_point<<-contador_direccion_point+1
  })
  
  ### MOVE POINT 
  observe({
    if (input$action_point=="move.event" & !is.null(data_of_click_remove_point$clickedMarkerRemove$lng)
        & !is.null(data_of_click_new_point$clickedPointNew$lng)){
      
      leafletProxy("mymapPointPatternEdition")%>%addCircleMarkers(lng=data_of_click_new_point$clickedPointNew$lng,
                                               lat=data_of_click_new_point$clickedPointNew$lat,
                                               radius=4,fillOpacity = 1,color="blue")
      
      print(paste("buscar_x_y"))
      print(head(Puntos_Actuales))
      print(data_of_click_remove_point$clickedMarkerRemove$lng)
      buscar_x=which((abs(Puntos_Actuales[,1]-data_of_click_remove_point$clickedMarkerRemove$lng)<10^(-10))==T)
      buscar_y=which((abs(Puntos_Actuales[,2]-data_of_click_remove_point$clickedMarkerRemove$lat)<10^(-10))==T)
      buscar_x_y=intersect(buscar_x,buscar_y)[1]

      if (!is.null(buscar_x_y)){
        if (!is.na(buscar_x_y)){
          Puntos_Mover<<-rbind(Puntos_Mover,c(buscar_x_y,data_of_click_new_point$clickedPointNew$lng,
                                           data_of_click_new_point$clickedPointNew$lat))
        }
      }

      data_of_click_remove_point$clickedMarkerRemove=NULL
      data_of_click_new_point$clickedPointNew=NULL
    }
  })
  
  ### REBUILDING THE POINT PATTERN
  
  ReBuild_PointPattern <- eventReactive(input$rebuildpoint, {
    
    ### UPDATE ZOOM LEVEL
    zoom_mapa_point<<-input$mymapPointPatternEdition_zoom
    bounds_mapa_point=input$mymapPointPatternEdition_bounds
    center_mapa_point<<-list(x=mean(c(bounds_mapa_point$east, bounds_mapa_point$west)),
                             y=mean(c(bounds_mapa_point$north, bounds_mapa_point$south)))
    
    for (i in c(1:nrow(Puntos_Mover))){
      Puntos_Actuales[Puntos_Mover[i,1],]<<-c(Puntos_Mover[i,2],Puntos_Mover[i,3])
    }
    
    aux_linnet=spatial_lines_lin
    
    print(input$utm_select)
    Puntos_Actuales_UTM=LONLAT2UTM(Puntos_Actuales,zone=input$utm_select)
    head(Puntos_Actuales_UTM)
    incidentes_ppp=ppp(x=Puntos_Actuales_UTM[,1],y=Puntos_Actuales_UTM[,2], 
                       c(object_global$domain$window$xrange[1]-1000,object_global$domain$window$xrange[2]+1000),
                       c(object_global$domain$window$yrange[1]-1000,object_global$domain$window$yrange[2]+1000))
    point_pattern_edited<<-lpp(incidentes_ppp,object_global$domain)
    marks(point_pattern_edited)<<-data.frame(Marcas_Actuales)
    ### CLEAN OBJECT
    Puntos_Mover<<-c()
    ### point_pattern_edited IS A lppx
    return(point_pattern_edited)
  })
  
})

