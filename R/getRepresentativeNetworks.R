#' Outputs representative networks for clades/subpopulations larger than a size filter (very small subpopulations are not considered in downstream analyses)
#' 
#' @param sampleIDs sampleID vector
#' @param SubpopSizeFilter the cutoff for small subpopulations. Smaller subpopulations have unstable covariance structure, so no network structure is calculated
#' @param num_networkEdge the number of edges to draw for each subpopulation mutual information network
#' @export


getRepresentativeNetworks<-function(sampleIDs, SubpopSizeFilter, num_networkEdge){
  
  newer_subpopulationLabels<-NULL
  inputMatrix_withSampleName<-NULL
  
  for(i in 1:length(sampleIDs)){
    #print(i)
    sampleID<-sampleIDs[i]
    load(paste0(sampleID,"_new_subpopulations_Representative_subpops.Rdata"))
    load(paste0(sampleID,"_dataMatrix.Rdata"))
    
    mainDir <- getwd()
    subDir <- paste0("/", sampleID, "_CladeNetworks")
    
    #output representative/clade networks
    dir.create(file.path(mainDir, subDir))
    setwd(paste0(mainDir, subDir))
    
    newlabels<-newer_subpopulationLabels  
    filteredClades<-names(table(newlabels))[table(newlabels)>SubpopSizeFilter]
    whetherToKeep<-newlabels %in% filteredClades
    newlabels<-newlabels[whetherToKeep]
    dataMatrix_filtered<-inputMatrix_withSampleName[whetherToKeep,]
    
    outputRepresentativeNetworks_topEdges(dataMatrix_filtered, newlabels, threshold=num_networkEdge)
    setwd(mainDir)
    
    save(dataMatrix_filtered, file=paste0(sampleID,"_dataMatrix_filtered.Rdata"))
    save(newlabels, file=paste0(sampleID,"_new_subpopulations_Representative_subpops_filtered.Rdata"))
    
  }  
  
}
