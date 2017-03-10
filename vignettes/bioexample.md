Biology Example
---------------------------

The PAC-MAN data analysis pipeline can be applied to mass-cytometry (CyTOF) data analysis. In this case, the user reads in the example data files (already saved as the Rdata format) subsetted from Bendall et al., 2011 and goes through the data analysis pipeline.

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
samplePass(sampleIDs, hyperrectangles=35, num_PACSupop=25, num_networkEdge=25, max.iter=50)
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
getRepresentativeNetworks(sampleIDs, SubpopSizeFilter=10, num_networkEdge=25)
```

Obtain annotations of subpopulations


```r
aggregateMatrix_withAnnotation<-annotateClades(sampleIDs, topHubs=4)
head(aggregateMatrix_withAnnotation)
```

```
##                        Annotation ClusterID SampleID   pPLCgamma2
## 1           pNFkB-pSTAT3-pH3-Ki67    clade1    Basal  0.868677090
## 2      IkBalpha-pP38-Ki67-pERK1.2    clade2    Basal  0.417650152
## 3      pNFkB-Ki67-pMAPKAPK2-pSHP2    clade3    Basal -0.009067766
## 4             Ki67-pP38-pS6-pSHP2   esc_1_1    Basal  0.477838191
## 5 pSTAT3-pNFkB-pERK1.2-pZAP70.Syk   esc_1_2    Basal  1.437341786
## 6    pSTAT5-pPLCgamma2-pNFkB-pP38   esc_1_3    Basal  1.053735073
##      pSTAT5       Ki67      pSHP2   pERK1.2 pMAPKAPK2 pZAP70.Syk    pSTAT3
## 1 1.4210433 1.60451708 1.04433647 1.4661011 1.5628676 1.43355321 2.1289258
## 2 0.7512015 0.75917867 0.60733563 0.9891308 0.9838982 0.68705268 1.4899395
## 3 0.1439336 0.59702795 0.06754460 0.1806816 0.2359121 0.05943879 0.1938544
## 4 1.0233990 2.76398013 0.52754587 0.6222940 0.9493720 0.62533761 0.9687189
## 5 2.0024983 6.45517278 1.46724050 1.6268664 4.4944990 1.81820683 3.3221809
## 6 6.6284754 0.03621155 0.03204475 1.1074446 0.1681345 0.17638086 0.5361319
##         pSLP     pNFkB  IkBalpha       pH3      pP38  pBtk.Itk       pS6
## 1 0.92427696 2.5118414 1.9303532 2.0931644 2.3369804 2.6054058 1.5696688
## 2 0.46984991 1.7826376 1.4261690 0.9924652 1.5957836 1.7199279 0.6538710
## 3 0.06218395 0.4461451 0.2065770 0.2719665 0.2130849 0.5528641 0.1526701
## 4 0.23311470 1.8449723 0.7085709 0.7539001 0.7906732 5.7284289 1.9941594
## 5 1.90283364 3.3438307 1.9513517 6.8553922 2.0497161 3.7951764 6.0227846
## 6 0.14897357 1.0628547 3.1909841 0.2512816 0.3274174 0.5179238 0.3013488
##      pSrcFK      pCREB      pCrkL count
## 1 2.9274172 1.80653721 1.03088589  1481
## 2 1.8988397 0.86030693 0.57582728   507
## 3 0.2939679 0.08732748 0.06956507   458
## 4 0.8047753 0.72321317 0.21448948   153
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
##                        Annotation ClusterID SampleID   pPLCgamma2
## 1           pNFkB-pSTAT3-pH3-Ki67    clade1    Basal  0.868677090
## 2      IkBalpha-pP38-Ki67-pERK1.2    clade2    Basal  0.417650152
## 3      pNFkB-Ki67-pMAPKAPK2-pSHP2    clade3    Basal -0.009067766
## 4             Ki67-pP38-pS6-pSHP2   esc_1_1    Basal  0.477838191
## 5 pSTAT3-pNFkB-pERK1.2-pZAP70.Syk   esc_1_2    Basal  1.437341786
## 6    pSTAT5-pPLCgamma2-pNFkB-pP38   esc_1_3    Basal  1.053735073
##      pSTAT5       Ki67      pSHP2   pERK1.2 pMAPKAPK2 pZAP70.Syk    pSTAT3
## 1 1.4210433 1.60451708 1.04433647 1.4661011 1.5628676 1.43355321 2.1289258
## 2 0.7512015 0.75917867 0.60733563 0.9891308 0.9838982 0.68705268 1.4899395
## 3 0.1439336 0.59702795 0.06754460 0.1806816 0.2359121 0.05943879 0.1938544
## 4 1.0233990 2.76398013 0.52754587 0.6222940 0.9493720 0.62533761 0.9687189
## 5 2.0024983 6.45517278 1.46724050 1.6268664 4.4944990 1.81820683 3.3221809
## 6 6.6284754 0.03621155 0.03204475 1.1074446 0.1681345 0.17638086 0.5361319
##         pSLP     pNFkB  IkBalpha       pH3      pP38  pBtk.Itk       pS6
## 1 0.92427696 2.5118414 1.9303532 2.0931644 2.3369804 2.6054058 1.5696688
## 2 0.46984991 1.7826376 1.4261690 0.9924652 1.5957836 1.7199279 0.6538710
## 3 0.06218395 0.4461451 0.2065770 0.2719665 0.2130849 0.5528641 0.1526701
## 4 0.23311470 1.8449723 0.7085709 0.7539001 0.7906732 5.7284289 1.9941594
## 5 1.90283364 3.3438307 1.9513517 6.8553922 2.0497161 3.7951764 6.0227846
## 6 0.14897357 1.0628547 3.1909841 0.2512816 0.3274174 0.5179238 0.3013488
##      pSrcFK      pCREB      pCrkL count subpop_proportion
## 1 2.9274172 1.80653721 1.03088589  1481           55.8900
## 2 1.8988397 0.86030693 0.57582728   507           19.1300
## 3 0.2939679 0.08732748 0.06956507   458           17.2800
## 4 0.8047753 0.72321317 0.21448948   153            5.7740
## 5 3.7489126 1.51700775 0.53621926    20            0.7547
## 6 0.5242967 0.63514127 0.38188762    31            1.1700
```




