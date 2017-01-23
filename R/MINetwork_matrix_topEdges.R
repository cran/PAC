#' @import parmigene
NULL

#' Mutual information network connection matrix generation (mrnet algorithm) using the parmigene package
#' 
#' @param dataMatrix data matrix with first column being the sample ID
#' @param threshold the number of edges to draw for each subpopulation mutual information network
#' @return the mutual information network connection matrix with top edges


MINetwork_matrix_topEdges<-function(dataMatrix, threshold){
  
  dataMatrix_mi<-knnmi.all(t(dataMatrix))
  dataMatrix_net<-mrnet(dataMatrix_mi)
  
  filterLevel<-dataMatrix_net[order(dataMatrix_net, decreasing=TRUE)][threshold]
  dataMatrix_net[dataMatrix_net<=filterLevel]<-0
  dataMatrix_net[dataMatrix_net>=filterLevel]<-1
  return(dataMatrix_net)
  
}