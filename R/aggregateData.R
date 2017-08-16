#' @importFrom dplyr %>% group_by summarise_each funs
NULL

#' Aggregates results from the clustering and merging step.
#' 
#' @param dataInput Data matrix, with first column being SampleID.
#' @param labelsInput cluster labels from PAC.
#' @return The aggregated data of \code{dataInput}, with average signal levels for all clusters and sample combinations.
#' @export


aggregateData<-function(dataInput, labelsInput){
  ClusterID<-NULL
  SampleID<-NULL
  
  data<-cbind(labelsInput, dataInput)
  colnames(data)[1]<-"ClusterID"
  colnames(data)[2]<-"SampleID"
  data<-data.frame(data)
  data$count<-1
  data_agg_intermediate<-data %>% group_by(ClusterID, SampleID) %>% summarise_each(funs(sum))
  data_agg_intermediate<-data.frame(data_agg_intermediate)
  data_agg<- data_agg_intermediate[order(data_agg_intermediate$SampleID),]
  data_agg[,-c(1,2,ncol(data_agg))]<-data_agg[,-c(1,2,ncol(data_agg))]/data_agg$count
  
  return(data_agg)
  
}