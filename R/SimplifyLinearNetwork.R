#' Simplifies a linear network object by merging the edges that meet at a second-degree vertex if some conditions on the angle they form and on their lengths are satisfied
#' 
#' This algorithm attempts to automatically reduce a linear network's complexity without altering its basic geometric configuration. The main objective of the algorithm is to merge the pairs of edges of the network that are connected by a second-degree vertex (with only two incident edges) into only one edge. Equivalently, this action means to join two vertex of the network whose path of connection only passes through another vertex of the network. 
#' 
#' @param network - A \code{linnet} object representing a linear network structure
#' @param Angle - An angle (in degrees, from 0 to 90) to adjust the level of simplification performed by the algorithm
#' @param Length - A length (in meters) to adjust the level of simplification performed by the algorithm
#' @param M - A k x 2 \code{matrix} that represents a combined effect for the \code{Angle} and \code{Length} parameters. First column represents a partition of the length space, whereas the second contains the angles considered in the simplification procedure when the lengths of both edges meeting in a second-degree vertex are below the threshold indicated in the first column
#' @return Returns a simplification of the linear network passed as an input (linnet class)
#' @examples 
#' library(SpNetPrep)
#' library(spatstat)
#' library(sp)
#' library(raster)
#' library(maptools)
#' network <- chicago$domain
#' Angle <- 10
#' Length <- 20
#' simplified_network_1 <- SimplifyLinearNetwork(network,Angle=Angle,Length=Length)
#' \dontrun{
#' M <- matrix(c(10,60,40,25),nrow=2)
#' simplified_network_2 <- SimplifyLinearNetwork(network,M=M)
#' }
#' @export
SimplifyLinearNetwork <- function(network,Angle=NULL,Length=NULL,M=NULL){
  edges_vertex=cbind(network$from,network$to)
  grado_2=which((vertexdegree(network)==2)==T)
  for (i in c(1:length(grado_2))){
    #print(i)
    ### A SECOND-DEGREE VERTEX, v, IS SELECTED
    ### EDGES THAT MEET AT THIS VERTEX ARE FOUND
    ### THE TWO OTHER VERTEX CONNECTED TO v ARE FOUND
    v=grado_2[i]
    buscar_1=which((edges_vertex[,1]==v)==T)
    buscar_2=which((edges_vertex[,2]==v)==T)
    ejes_eliminar=c(buscar_1,buscar_2)
    puntos_unir=c(edges_vertex[buscar_1,],edges_vertex[buscar_2,])
    puntos_unir=puntos_unir[puntos_unir!=v]
    
    ### ANGLES THAT FORM THE TWO EDGES THAT MEET AT v
    
    resta1=c(network$vertices$x[puntos_unir[1]],network$vertices$y[puntos_unir[1]])-
      c(network$vertices$x[v],network$vertices$y[v])
    angulo1=atan(resta1[2]/resta1[1])*(180/pi)
    
    resta2=c(network$vertices$x[puntos_unir[2]],network$vertices$y[puntos_unir[2]])-
      c(network$vertices$x[v],network$vertices$y[v])
    angulo2=atan(resta2[2]/resta2[1])*(180/pi)
    
    ### LENGTHS OF THE TWO EDGES THAT MEET AT v
    
    longitud1=pointDistance(c(network$vertices$x[v],network$vertices$y[v]),
                            c(network$vertices$x[puntos_unir[1]],network$vertices$y[puntos_unir[1]]),
                            lonlat = F)
    longitud2=pointDistance(c(network$vertices$x[v],network$vertices$y[v]),
                            c(network$vertices$x[puntos_unir[2]],network$vertices$y[puntos_unir[2]]),
                            lonlat = F)
    
    ### MODIFICATION DEPENDING ON THE INPUT INTRODUCED FOR CONTROLLING THE ACTION OF THE 
    ### ALGORITHM. FIRST CONDITION REPRESENTS THE SIMPLEST CASE, WHEN M=NULL. SECOND CONDITION 
    ### IS FOR A CHOICE OF M
    if (is.null(M)){
      condicion_angulo=abs(angulo1-angulo2)<=Angle | (180-abs(angulo1-angulo2))<=Angle
      condicion_longitud=longitud1<=Length | longitud2<=Length
      ### IF BOTH CONDITIONS HOLD THE EDGES ARE MERGED AND THE SECOND-DEGREE VERTEX REMOVED
      if (condicion_angulo & condicion_longitud){
        edges_vertex=edges_vertex[-ejes_eliminar,]
        edges_vertex=rbind(edges_vertex,c(puntos_unir[1],puntos_unir[2]))
      }
    } else{
      comprobar_l1=longitud1<=M[,1]
      comprobar_l2=longitud2<=M[,1]
      comprobar_producto=comprobar_l1*comprobar_l2
      condicion=which((comprobar_producto)==T)[1]
      ### IF condicion IS NA THEN BOTHE EDGES VIOLATE LENGTH CONDITION, NO MORE CHECKS ARE REQUIRED
      if (!is.na(condicion)){
        Angle_aux=M[condicion,2]
        condicion_angulo=abs(angulo1-angulo2)<=Angle_aux | (180-abs(angulo1-angulo2))<=Angle_aux
        ### IF CONDITION ON THE ANGLE HOLDS THE EDGES ARE MERGED AND THE SECOND-DEGREE VERTEX REMOVED
        if (condicion_angulo){
          edges_vertex=edges_vertex[-ejes_eliminar,]
          edges_vertex=rbind(edges_vertex,c(puntos_unir[1],puntos_unir[2]))
        }
      }
    }
  }
  ### NEW NETWORK IS CONSTRUCTED
  indice=0
  lineas=list()
  for (i in c(1:nrow(edges_vertex))){
    indice=indice+1
    linea=Lines(list(Line(rbind(c(network$vertices$x[edges_vertex[i,1]],
                                  network$vertices$y[edges_vertex[i,1]]),
                                c(network$vertices$x[edges_vertex[i,2]],
                                  network$vertices$y[edges_vertex[i,2]])))),
                toString(indice))
    lineas[[indice]]=linea
  }
  lineas <- SpatialLines(lineas)
  lineas <- as.linnet.SpatialLines(lineas)
  return(lineas)
}

