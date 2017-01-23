#' @import parmigene
#' @importFrom igraph graph.adjacency simplify
NULL

#' Plots mutual information network (using mrnet algorithm) connection using the parmigene package
#' @param dataMatrix data matrix with first column being the sample ID
#' @param threshold the number of edges to draw for each subpopulation mutual information network

MINetworkPlot_topEdges<-function(dataMatrix, threshold){
  igraph_opt<-NULL
  layout.circle<-NULL
  
  dataMatrix_mi<-knnmi.all(t(dataMatrix))
  dataMatrix_net<-mrnet(dataMatrix_mi)
  
  filterLevel<-dataMatrix_net[order(dataMatrix_net, decreasing=TRUE)][threshold]
  dataMatrix_net[dataMatrix_net<filterLevel]<-0
  
  g_dataMatrix_net <- graph.adjacency(dataMatrix_net, weighted=TRUE, mode="undirected")
  g_dataMatrix_net<-simplify(g_dataMatrix_net, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb=igraph_opt("edge.attr.comb"))
  networkPlot<-plot(g_dataMatrix_net, layout=layout.circle(g_dataMatrix_net))
}
