
#Function for finding if a splice variant has background-specific expression:

PSI = matrix of PSI values for splice variants (row names), one column per sample (column names).

````r

#                                       HA351CR1R1  HA351CR1R2  HA351CR1R3  HA351CR1R4
#HanXRQChr01g0012571_2_3_IR_na:na:+:IR  1.0000000           1         0.9           1
#HanXRQChr10g0310601_1_2_IR_na:na:+:IR  0.1071429 0.111111111           0 0.057142857
#HanXRQChr10g0310601_7_8_IR_na:na:+:IR  0.1333333 0.099099099 0.082644628 0.142857143
#HanXRQChr10g0285131_3_4_IR_na:na:+:IR  0.3142857 0.269230769 0.363636364 0.461538462
#HanXRQChr10g0285131_5_6_IR_na:na:+:IR  1.0000000           1 0.526315789           1
#HanXRQChr11g0322961_3_4_IR_na:na:+:IR  0.0000000           0         0.5         0.1

#SAM208 = matrix of SAM208 sample PSI values for the age~tissue combination of interest
#SAM004 = matrix of SAM004 sample PSI values for the age~tissue combination of interest

#SAM004 and #SAM208 are in the same order and have the same row names!
find_DAS_significance <- function(SAM208, SAM004){
  
  #find PSI means for each splice variant in both backgrounds:
  SAM208[] <- sapply(SAM208, as.numeric)
  SAM208$avg <- rowMeans(SAM208)
  SAM004[] <- sapply(SAM004, as.numeric)
  SAM004$avg <- rowMeans(SAM004)
  
  #for each splice variant:
  for(row in 1:nrow(SAM208)){
    splice <- rownames(SAM208[row,])
    
    #make a PSI vector, not including the mean, for the splice variant of interestin both genetic backgrounds:
    v208 <- as.numeric(as.character(SAM208[row,]))
    v208 <- v208[-5]
    
    v004 <- as.numeric(as.character(SAM004[row,]))
    v004 <- v004[-5]
    
    if(is.numeric(v208) & is.numeric(v004)){ tryCatch({
      t<- (t.test(v208, v004, alternative = "two.sided", conf.level = 0.95))
      p <- t$p.value
    }, error=function(e){p <- 1})} else { p <- NA}
    
    out <- as.data.frame(t(data.frame(c(splice, p))))
    
    rownames(out) <- c()
    #tagP <- paste("p_",age,sep="")
    colnames(out) <- c("Splice_Variant", "p")
    
    if(exists("p_val_out") == F){
      p_val_out <- out
    } else (p_val_out <- rbind(p_val_out,out))
    
  }
  
  p_val_out$FDR <- p.adjust(p_val_out$p, method = "fdr")
  return(p_val_out)
}

#Function returns a two-column df of splice variants and their p-value + p.Adj for background specificity. 


````

