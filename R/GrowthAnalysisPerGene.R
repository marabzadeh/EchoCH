#' GrowthAnalysisPerGene
#'
#' Growth analysis per-gene per-treatment
#'
#' @param PerVarGrowth.df a data.frame which is the output of GrowthAnalysisPerVariant
#' @param time_threshold a threshhold on the timepoints (in days)
#' @return   PerGeneGrowth.df : dataframe reporting the VAF-changes-per-month for each gene per treatment and other statistics
#' @examples
#' PerGeneGrowth.df <- GrowthAnalysisPerGene(PerVarGrowth.df, time_threshold)

GrowthAnalysisPerGene <- function(PerVarGrowth.df, time_threshold){

  file.df <- PerVarGrowth.df

  file.df$VAFChangePerMonth [ which(is.na(file.df$VAFChangePerMonth) == T) ] <- 0

  ran <- dim(table(file.df$GENE, file.df$treatment) )
  tab.ls <- table(file.df$GENE, file.df$treatment)


  gene.ls <- vector (mode = "character", length = ran[1]*ran[2]  )
  trt.ls <- vector (mode = "character", length = ran[1]*ran[2]  )
  mean.ls <- vector (mode = "integer", length = ran[1]*ran[2]  )
  median.ls <- vector (mode = "integer", length = ran[1]*ran[2]  )
  min.ls <-   vector (mode = "integer", length = ran[1]*ran[2]  )
  max.ls <-   vector (mode = "integer", length = ran[1]*ran[2]  )
  SD.ls <-   vector (mode = "integer", length = ran[1]*ran[2]  )
  SE.ls <-   vector (mode = "integer", length = ran[1]*ran[2]  )
  mutNum.ls <-vector (mode = "integer", length = ran[1]*ran[2]  )
  ptNum.ls <-vector (mode = "integer", length = ran[1]*ran[2]  )
  meanSE <-vector (mode = "integer", length = ran[1]*ran[2]  )
  datapoints <-vector (mode = "list", length = ran[1]*ran[2]  )

  count <- 1

  for (i in 1: ran[1]){
    for (j in 1:ran[2]){

      file.df.sub <- file.df[file.df$status_1 == "pre", ]
      file.df.sub <- file.df.sub[as.numeric(file.df.sub$DeltaTime) <= time_threshold, ]

      file.df.sub <- file.df.sub[file.df.sub$GENE == row.names(tab.ls)[i], ]
      file.df.sub <- file.df.sub[file.df.sub$treatment == colnames(tab.ls)[j], ]
      #file.df.sub$dVAF <- as.numeric(file.df.sub$finalAF_2) - as.numeric(file.df.sub$finalAF_1)

      if (dim (file.df.sub)[1] !=0 ){

        id.lst <- unique(file.df.sub$id )
        ls.max <- NULL
        ls.mut <- NULL

        if ( length(id.lst) == 1){
          file.df.sub.sub <- file.df.sub[ file.df.sub$id == id.lst, ]
          #ls.max <- max(file.df.sub.sub$dVAF)
          ls.max <- max(file.df.sub.sub$VAFChangePerMonth)
          ls.mut <- sort(table(file.df.sub.sub$order_2) , decreasing = T)[1]
        }else{
          for (k in 1:length(id.lst)){

            file.df.sub.sub <- file.df.sub[ file.df.sub$id == id.lst[k], ]
            #ls.max <- append (ls.max, max(file.df.sub.sub$dVAF) )
            ls.max <- append (ls.max, max(file.df.sub.sub$VAFChangePerMonth) )

            ls.mut <- append (ls.mut, sort(table(file.df.sub.sub$order_2) , decreasing = T)[1] )
          }
        }#end else > 1 pt

        ptNum.ls[count] <- length(id.lst)
        mutNum.ls [count] <- sum(ls.mut)

        gene.ls[count] <- row.names(tab.ls)[i]
        trt.ls[count] <- colnames(tab.ls)[j]

        #CalcStatRes <-  CalcStat(ls.max)

        datapoints[[count]] <- paste(ls.max,collapse = ",")
        mean.ls[count] <- mean(ls.max)
        median.ls[count] <- median(ls.max)

        min.ls[count] <- min(ls.max)
        max.ls[count] <- max(ls.max)

        SD.ls[count] <- sd(ls.max)
        SE.ls[count] <- SD.ls[count] / sqrt( mutNum.ls[count] )
        meanSE[count] <- mean.ls[count]- SE.ls[count]

        #if ( row.names(tab.ls)[i] == "ASXL1"){
        #  print(count)
        #  print(i)
        #  print(j)}

        count <- count+1

      }## end if

    }## end for j
  }## end for i

  datapoints <- unlist(datapoints)
  nn <- length(datapoints)

  res.df <- data.frame(gene.ls[1:nn], trt.ls[1:nn], mean.ls[1:nn],
                       median.ls[1:nn],
                       min.ls[1:nn], max.ls[1:nn],
                       SD.ls[1:nn], SE.ls[1:nn], ptNum.ls[1:nn], mutNum.ls[1:nn], meanSE[1:nn], datapoints)
  res.df <- res.df[res.df$gene.ls != "",]

  res.df <- res.df[order(res.df$meanSE, decreasing = T), ]

  colnames (res.df) <- c("gene",	"treatment",	"mean_CHvaf_change_per_month",	"median_CHvaf_change_per_month",
                         "min_CHvaf_change_per_month"	, "max_CHvaf_change_per_month",
                         "standard_deviation",	"standard_error",	"number_patients",	"number_mutations",
                         "meanSE", "datapoints")

  PerGeneGrowth.df <- res.df

  return(PerGeneGrowth.df)

}
