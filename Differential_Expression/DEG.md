# R script for finding DEG
**Packages: EdgeR, Limma, tidyverse, reshape2**



### Files Needed: 
1. tissue-specific meta data (here is leaf only):

````
#         Background Tissue Treatment Age Ann_01_01 Ann_05_01
#208CL1R1     SAM208   Leaf   Control  16         0         2
#208CL1R2     SAM208   Leaf   Control  16         0         2
#208CL1R3     SAM208   Leaf   Control  16         0         2
#208CL1R4     SAM208   Leaf   Control  16         0         2
#208CL2R1     SAM208   Leaf   Control  24         0         2
#208CL2R2     SAM208   Leaf   Control  24         0         2
````

2. raw counts:

````
#                 208CL1R1 208CL1R2 208CL1R3 208CL1R4 208CL2R1 208CL2R2 208CL2R3 208CL2R4 208CL3R1 208CL3R2
#HanXRQCPg0579731      792      683      867      554      347      286      341      447      417      874
#HanXRQCPg0579741       60       98      114       50       30      157       60      218       27       27
#HanXRQCPg0579761      838     1159     1243      971     1340     4686     1006     8456      841     1095
#HanXRQCPg0579801      125      140      188      155       74      254       94      281       28       48
#HanXRQCPg0579811        1        0        2        0        0        0        0        1        0        5
#HanXRQCPg0579821        5        6        6        0        3        5        3        3        3       12
#                 208CL3R3 208CL3R4 4CL1R1 4CL1R2 4CL1R3 4CL1R4 4CL2R1 4CL2R2 4CL2R3 4CL2R4 4CL3R1 4CL3R2 4CL3R3
#HanXRQCPg0579731      393      862    602    365    921    415    276    341    466    312    706    608    501
#HanXRQCPg0579741      103       15     82     70     78     65     58     79    136     91     46     17     34
#HanXRQCPg0579761     2332      763   1900   1722   1619   2160   1408   2431   2367   2591   1378    698   1141
#HanXRQCPg0579801      153       28    131    130     93    166     93    163    235    190     93     36     81
#HanXRQCPg0579811        0        1      0      0      3      2      1      1      0      0      1      2      1
#HanXRQCPg0579821        4        7      7      1     16      7      6      2      6      7      0      4      6
#                 4CL3R4
#HanXRQCPg0579731    537
#HanXRQCPg0579741     15
#HanXRQCPg0579761   1108
#HanXRQCPg0579801     56
#HanXRQCPg0579811      0
#HanXRQCPg0579821      4
````

3. Choose a cpm threshold (cpm below which genes should be filtered out)

````
cpm_thresh <- 1
````

4. Specify tissue

````
tissue <- "leaf"
````

#### Find DEG by tissue ~ age combination between the two SAM lines using EdgeR/Limma:

````
All_Genes_Limma <- function(meta, counts, cpm_thresh, tissue){
  
  meta <- meta
  
  for(age in unique(meta$Age)){
    tissue <- tissue
    temp_meta <- meta[meta$Age == age,]
    #rownames(temp_meta) <- c()
    #temp_meta <- temp_meta %>% column_to_rownames("Individual")
    temp_counts <- counts[, colnames(counts)%in%rownames(temp_meta)]
    
    tCPM <- cpm(temp_counts)
    thresh <- tCPM > cpm_thresh
    keep <- rowSums(thresh) >= 4
    counts.keep <- temp_counts[keep,]
    y <- DGEList(counts.keep, group = temp_meta$Background)
    y <- calcNormFactors(y)
    
    #calculate dispersion:
    dgeObj <- estimateCommonDisp(y)
    dgeObj <- estimateGLMTrendedDisp(dgeObj)
    dgeObj <- estimateTagwiseDisp(dgeObj)
    
    #Make deisgn_Matrix
    f <- factor(temp_meta$Background, levels = c("SAM004", "SAM208"))
    designMat <- model.matrix(~0+f)
    colnames(designMat) <- c("SAM004", "SAM208")
    rownames(designMat) <- rownames(temp_meta)
    
    #genewise glms method:
    
    #Find differentially expressed genes by inversion:
    fit <- glmFit(dgeObj, designMat)
    names(fit)
    head(coef(fit))
    PvsV <- makeContrasts(SAM004-SAM208, levels=designMat)
    lrt.pVsV <- glmLRT(fit, contrast=PvsV)
    genes <- as.data.frame(topTags(lrt.pVsV, n = nrow(lrt.pVsV), adjust.method="BH"))
    
    genes <- genes[c(1,5)]
    genes <- rownames_to_column(genes, "Gene")
    
    colnames(genes)[2:3] <- c(paste("d", age, "_LogFC",sep=""), paste("d",age,"_FDR",sep=""))
    
    if(exists("Out") == F){
      Out <- genes
    } else {Out <- merge(Out, genes, by ="Gene", all = T)}
  }
  
  return(Out)
}
````

**Function Outputs the logFC and FDR-corrected p-value for each gene in leaves across the three ages (16, 24, 30 days)**

````
#                      Gene  d16_LogFC      d16_FDR   d24_LogFC     d24_FDR  d30_LogFC    d30_FDR
#1 HanXRQChr00c0004g0571011 -0.8018670 6.719667e-02 -1.01772615 0.002691688 -0.8103942 0.00535298
#2 HanXRQChr00c0016g0571081  0.6830112 8.514026e-02 -0.16637376 0.724782872  0.4934317 0.41412832
#3 HanXRQChr00c0020g0571101  0.4332272 1.994385e-01  1.01116961 0.093119749  0.4893597 0.14714535
#4 HanXRQChr00c0028g0571141  0.0510072 8.162738e-01 -0.01561328 0.957830705 -0.2729970 0.22902954
#5 HanXRQChr00c0029g0571181 -1.1031446 2.101241e-07  0.28108343 0.689375452  0.2195302 0.67283551
#6 HanXRQChr00c0036g0571221 -0.3953910 5.642824e-02 -0.64421789 0.042856452 -0.6886883 0.01415004
````

**DEG enrichment within an inversion can be checked with the following function:**

````
significance <- function(All_DEG, inverted_genes, all_genes){
  #All_DEG = character vector for gene names of all DEG
  #inverted_genes = character vector for gene names of all inverted genes
  #all_genes = numerical all genes in the genome. IN this case, nrow(xrq)
  
  inv_DEG <- length(All_DEG[All_DEG%in%inverted_genes])   #Inverted DEG
  noninv_DEG <- length(All_DEG)-inv_DEG
  
  inv_nonDEG <- length(inverted_genes)-inv_DEG
  noninv_nonDEG <-  all_genes - noninv_DEG
  
  #Two-sided Fisher's exact test:
  #_______________#______________#___________________#
  #               #    DfExpr    #    not DfExpr     #
  #_______________#______________#___________________#
  #   Inverted   #1       ?      #3        ?         #
  #_______________#______________#___________________#
  # non-Inverted  #2      ?      #4        ?         #
  #_______________#______________#___________________#
  #
  
  fisher_res <- stats::fisher.test(matrix(c((inv_DEG),(noninv_DEG),(inv_nonDEG), (noninv_nonDEG)),nrow=2), y = NULL, alternative = "greater", conf.level = 0.95)
  #p-value for Ann_01_01 is 0.75 
  fisher_p <- fisher_res$p.value
  
  return(fisher_p)
}
````

The function checks with a Fisher's Exact Test if the number of inverted genes is greater than what's expected under random chance. Checks if [inverted DEG / non-DEG inverted genes] is significantly greater than [non-inverted DEG / non-DEG non-inverted genes]










