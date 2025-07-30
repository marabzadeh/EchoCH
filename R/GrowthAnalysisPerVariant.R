#' GrowthAnalysisPerVariant
#'
#' Growth analysis per-variant per-patient for each time-point
#'
#' @param CH.df a data.frame with CH mutations with the following columns::
#'                                              patient_id (combination of number and letter), gene, treatment, amino_acid_change,
#'                                              depth, vaf, temporal_order, pre/post_treatment,
#'                                             days_from_treatment_start_to_sample_collection )
#' @return PerVarGrowth.df : dataframe reporting the VAF-changes-per-month and other statistics
#' @examples
#' PerVarGrowth.df <- GrowthAnalysisPerVariant(CH.df)
GrowthAnalysisPerVariant <- function(CH.df){

  CH.df$depth <- as.numeric(CH.df$depth)
  CH.df$SampVar <- paste0(CH.df$patient_id, "-" , CH.df$amino_acid_change)

  uniqueLoop <- (unique(CH.df$SampVar))

  # id - treatment - gene
  # order status AF DP
  # order status AF DP
  # timeDef

  uniqueLoop <- (unique(CH.df$SampVar))

  df.new <- data.frame(id=c() , treatment=c(), GENE= c(),
                       order_1=c(), status_1=c(), finalAF_1=c(), totalDepth_1=c(),
                       order_2=c(), status_2=c(), finalAF_2=c(), totalDepth_2=c(),
                       DeltaTime=c())

  for (i in 1: length(uniqueLoop) ){

    tmp.df <- CH.df[CH.df$SampVar == uniqueLoop[i], ]

    if (dim(tmp.df)[1] == 2 ){

      firstInx <- which(tmp.df$temporal_order == 1 )
      SecInx <- which(tmp.df$temporal_order == 2 )

      tmp.mk <- data.frame(tmp.df$patient_id[firstInx], tmp.df$treatment[firstInx], tmp.df$gene[firstInx],
                           tmp.df$temporal_order[firstInx], tmp.df$`pre/post_treatment`[firstInx], tmp.df$vaf[firstInx], tmp.df$depth[firstInx],
                           tmp.df$temporal_order[SecInx], tmp.df$`pre/post_treatment`[SecInx], tmp.df$vaf[SecInx], tmp.df$depth[SecInx],
                           as.numeric(tmp.df$days_from_treatment_start_to_sample_collection[SecInx] ) - as.numeric(tmp.df$days_from_treatment_start_to_sample_collection[firstInx]) )
      colnames(tmp.mk) <- c("id" , "treatment", "GENE",
                            "order_1", "status_1", "finalAF_1", "totalDepth_1",
                            "order_2", "status_2", "finalAF_2", "totalDepth_2",
                            "DeltaTime")
      df.new <- rbind(df.new, tmp.mk)
    }else{

      for(j in 1:(max(tmp.df$temporal_order)-1) ) {
        for (k in (j+1):max(tmp.df$temporal_order)){

          firstInx <- which(tmp.df$temporal_order == j )
          SecInx <- which(tmp.df$temporal_order == k )

          tmp.mk <- data.frame(tmp.df$patient_id[firstInx], tmp.df$treatment[firstInx], tmp.df$gene[firstInx],
                               tmp.df$temporal_order[firstInx], tmp.df$`pre/post_treatment`[firstInx], tmp.df$vaf[firstInx], tmp.df$depth[firstInx],
                               tmp.df$temporal_order[SecInx], tmp.df$`pre/post_treatment`[SecInx], tmp.df$vaf[SecInx], tmp.df$depth[SecInx],
                               as.numeric(tmp.df$days_from_treatment_start_to_sample_collection[SecInx] ) - as.numeric(tmp.df$days_from_treatment_start_to_sample_collection[firstInx]) )
          colnames(tmp.mk) <- c("id" , "treatment", "GENE",
                                "order_1", "status_1", "finalAF_1", "totalDepth_1",
                                "order_2", "status_2", "finalAF_2", "totalDepth_2",
                                "DeltaTime")
          df.new <- rbind(df.new, tmp.mk)
        }
      }
    }
  }

  #(G*H+L*M)/(H+M)
  df.new$pHat <- ( as.numeric(df.new$finalAF_1)*as.numeric(df.new$totalDepth_1) +  as.numeric(df.new$finalAF_2)*as.numeric(df.new$totalDepth_2) )/( as.numeric(df.new$totalDepth_1) + as.numeric(df.new$totalDepth_2) )
  #SQRT(N*(1-N)*(1/H+1/M))
  df.new$sigma <- sqrt(df.new$pHat*(1-df.new$pHat)*(1/as.numeric(df.new$totalDepth_1) + 1/as.numeric(df.new$totalDepth_2) )   )
  #df.new$sigma2 <- sqrt( (df.new$pHat*(1-df.new$pHat) ) / ( as.numeric(df.new$totalDepth_1) + as.numeric(df.new$totalDepth_2) )   )
  #df.new$sigma2 - df.new$sigma

  #L-G
  df.new$dVAF <- as.numeric(df.new$finalAF_2) - as.numeric(df.new$finalAF_1)
  #(L-G)/SQRT(N*(1-N)*(1/H+1/M))
  df.new$z <- df.new$dVAF / df.new$sigma

  df.new$DeltaTimeMonth <- df.new$DeltaTime /30
  #Q/T*30
  df.new$SigmaChangePerMonth <- df.new$z / df.new$DeltaTimeMonth
  df.new$VAFChangePerMonth <- df.new$z / df.new$DeltaTimeMonth * df.new$sigma

  PerVarGrowth.df <- df.new

  return (PerVarGrowth.df)

}
