#' Point pattern on a road network simulating traffic accidents
#'
#' An \code{lpp} object generated with a Poisson process of intensity 0.01 over the \code{SampleNetwork} road structure. The pattern is marked with three categorical variables: \code{accident.date}, \code{collision} and \code{vehicles}, which were also randomly constructed
#'
"SamplePointPattern"
#' Linear network representing a road structure
#'
#' An \code{SpatialLines} representing a road structure 
#'
"SampleNetwork"
#' Linear network representing a directed road structure
#'
#' An \code{SpatialLinesDataFrame} representing a directed road structure. The \code{data.frame} attached has three columns, \code{V1}, \code{V2} and \code{Dir} that determine the direction of the network. The column \code{Dir} indicates the type of flow existing between \code{V1} and \code{V2}: none (0), from \code{V1} to \code{V2} (1), from \code{V2} to \code{V1} (-1) or bidirectional (2).  
#'
"SampleDirectedNetwork"