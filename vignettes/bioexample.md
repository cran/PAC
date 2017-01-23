Biology Example
---------------------------

The PAC-MAN data analysis pipeline can be applied to mass-cytometry (CyTOF) data analysis. In this case, the user reads in the example data files (already saved as the Rdata format) from Bendall et al., 2011 and goes through the data analysis pipeline.

Load the required R packages


```r
library(NMF) #for plotting heatmap; the user can choose to use other heatmap packages as well

library(PAC)
```

Construct the sampleIDs vector to analyze the data


```r
sampleIDs<-c("Basal", "BCR", "IL7")
```

Partition, cluster into desired number of subpopulations, and output subpopulation mutual information networks


```r
samplePass(sampleIDs, hyperrectangles=35, num_PACSupop=25, num_networkEdge=50, max.iter=50)
```

```
## Input Data: 2650 by 18
## Partition method: Discrepancy based partition
## Maximum level: 35
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

```
## Input Data: 3537 by 18
## Partition method: Discrepancy based partition
## Maximum level: 35
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

```
## Input Data: 3813 by 18
## Partition method: Discrepancy based partition
## Maximum level: 35
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

Multiple Alignments of Networks


```r
clades_network_only<-MAN(sampleIDs, num_PACSupop=25, smallSubpopCutoff=100, k_clades=5)
```

Refine the PAC labels with multiple alignments of networks representative labels for clades


```r
refineSubpopulationLabels(sampleIDs,clades_network_only, expressionGroupClamp=5)
```

Draw clade/representative mutual information networks


```r
getRepresentativeNetworks(sampleIDs, SubpopSizeFilter=10, num_networkEdge=50)
```

Obtain annotations of subpopulations


```r
aggregateMatrix_withAnnotation<-annotateClades(sampleIDs, topHubs=4)
head(aggregateMatrix_withAnnotation)
```

```
##                     Annotation ClusterID SampleID   pPLCgamma2    pSTAT5
## 1      pNFkB-pSTAT3-pSrcFK-pH3    clade1    Basal  0.895183514 1.4577365
## 2  IkBalpha-pERK1.2-pNFkB-pP38    clade2    Basal  0.406696788 0.7133375
## 3     pNFkB-pH3-pMAPKAPK2-pP38    clade3    Basal -0.009067766 0.1439336
## 4 pSrcFK-IkBalpha-pSTAT3-pNFkB   esc_1_1    Basal  0.710704942 1.4471062
## 5      pERK1.2-pSTAT3-pP38-pH3   esc_1_2    Basal  1.437341786 2.0024983
## 6     pNFkB-pSTAT5-pSTAT3-pP38   esc_1_3    Basal  1.053735073 6.6284754
##         Ki67      pSHP2   pERK1.2 pMAPKAPK2 pZAP70.Syk    pSTAT3
## 1 1.64048700 1.09366726 1.5035418 1.5821762 1.51613641 2.1894451
## 2 0.60229621 0.51935971 0.9303466 0.9590697 0.60039036 1.3581764
## 3 0.59702795 0.06754460 0.1806816 0.2359121 0.05943879 0.1938544
## 4 3.57887762 0.87461857 1.0226867 1.4136084 0.98768716 1.6250415
## 5 6.45517278 1.46724050 1.6268664 4.4944990 1.81820683 3.3221809
## 6 0.03621155 0.03204475 1.1074446 0.1681345 0.17638086 0.5361319
##         pSLP     pNFkB IkBalpha       pH3      pP38  pBtk.Itk       pS6
## 1 0.97219644 2.5943881 1.916078 2.2105594 2.3501491 2.5445134 1.6190865
## 2 0.39876320 1.6165028 1.463709 0.8577937 1.5937568 1.8541578 0.6370273
## 3 0.06218395 0.4461451 0.206577 0.2719665 0.2130849 0.5528641 0.1526701
## 4 0.56318931 2.4699900 1.214800 1.3646863 1.4552160 6.0029834 2.3336508
## 5 1.90283364 3.3438307 1.951352 6.8553922 2.0497161 3.7951764 6.0227846
## 6 0.14897357 1.0628547 3.190984 0.2512816 0.3274174 0.5179238 0.3013488
##      pSrcFK      pCREB      pCrkL count
## 1 3.0377862 1.90677482 1.08460926  1326
## 2 1.7706313 0.76626335 0.49825585   641
## 3 0.2939679 0.08732748 0.06956507   458
## 4 1.4842993 1.16523341 0.53982438   174
## 5 3.7489126 1.51700775 0.53621926    20
## 6 0.5242967 0.63514127 0.38188762    31
```


Obtain heatmap input and plot heatmap


```r
cladeProportionMatrix<-heatmapInput(aggregateMatrix_withAnnotation)
aheatmap(as.matrix(cladeProportionMatrix))
```

![plot of chunk Heatmap input and plot heatmap](figure/Heatmap input and plot heatmap-1.png)


Append subpopulation proportions for each sample in the annotation matrix


```r
annotationMatrix_prop<-annotationMatrix_withSubpopProp(aggregateMatrix_withAnnotation)
head(annotationMatrix_prop)
```

```
##                     Annotation ClusterID SampleID   pPLCgamma2    pSTAT5
## 1      pNFkB-pSTAT3-pSrcFK-pH3    clade1    Basal  0.895183514 1.4577365
## 2  IkBalpha-pERK1.2-pNFkB-pP38    clade2    Basal  0.406696788 0.7133375
## 3     pNFkB-pH3-pMAPKAPK2-pP38    clade3    Basal -0.009067766 0.1439336
## 4 pSrcFK-IkBalpha-pSTAT3-pNFkB   esc_1_1    Basal  0.710704942 1.4471062
## 5      pERK1.2-pSTAT3-pP38-pH3   esc_1_2    Basal  1.437341786 2.0024983
## 6     pNFkB-pSTAT5-pSTAT3-pP38   esc_1_3    Basal  1.053735073 6.6284754
##         Ki67      pSHP2   pERK1.2 pMAPKAPK2 pZAP70.Syk    pSTAT3
## 1 1.64048700 1.09366726 1.5035418 1.5821762 1.51613641 2.1894451
## 2 0.60229621 0.51935971 0.9303466 0.9590697 0.60039036 1.3581764
## 3 0.59702795 0.06754460 0.1806816 0.2359121 0.05943879 0.1938544
## 4 3.57887762 0.87461857 1.0226867 1.4136084 0.98768716 1.6250415
## 5 6.45517278 1.46724050 1.6268664 4.4944990 1.81820683 3.3221809
## 6 0.03621155 0.03204475 1.1074446 0.1681345 0.17638086 0.5361319
##         pSLP     pNFkB IkBalpha       pH3      pP38  pBtk.Itk       pS6
## 1 0.97219644 2.5943881 1.916078 2.2105594 2.3501491 2.5445134 1.6190865
## 2 0.39876320 1.6165028 1.463709 0.8577937 1.5937568 1.8541578 0.6370273
## 3 0.06218395 0.4461451 0.206577 0.2719665 0.2130849 0.5528641 0.1526701
## 4 0.56318931 2.4699900 1.214800 1.3646863 1.4552160 6.0029834 2.3336508
## 5 1.90283364 3.3438307 1.951352 6.8553922 2.0497161 3.7951764 6.0227846
## 6 0.14897357 1.0628547 3.190984 0.2512816 0.3274174 0.5179238 0.3013488
##      pSrcFK      pCREB      pCrkL count subpop_proportion
## 1 3.0377862 1.90677482 1.08460926  1326           50.0400
## 2 1.7706313 0.76626335 0.49825585   641           24.1900
## 3 0.2939679 0.08732748 0.06956507   458           17.2800
## 4 1.4842993 1.16523341 0.53982438   174            6.5660
## 5 3.7489126 1.51700775 0.53621926    20            0.7547
## 6 0.5242967 0.63514127 0.38188762    31            1.1700
```




