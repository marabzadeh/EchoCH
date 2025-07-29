#' FitnessAnalysis
#'
#' Calculate Fitness likelihood for each variant in each patient
#'
#' @param CH.df a data.frame with CH mutations with the following columns::
#'                                              patient_id, gene, treatment, amino_acid_change,
#'                                              depth, vaf, temporal_order, pre/post_treatment,
#'                                             days_from_treatment_start_to_sample_collection )
#' @param Neff_Pr_List Neff_Pr_List[[1]]: a dataframe containing Neff results for each patient;
#'                     Neff_Pr_List[[2]]: a dataframe containing Pr results for each patient
#'
#' @param SigOrderOfMagnitude_threshold a threshhold on the order of magnitude of the likelihood of the mutation fitness
#'                                        (in comparison to the lowest likelihood in that patient)
#' @return prList.df : dataframe which reports the fitness of each mutation in each patient and the dynamic (pos/neg/noSelection)
#' @examples
#' prList.df <- FitnessAnalysis(CH.df, Neff_Pr_List)

FitnessAnalysis <- function (CH.df, Neff_Pr_List, SigOrderOfMagnitude_threshold = 5  ){

  Pr.listS <- Neff_Pr_List[[2]]

  CH.df.tmp <- CH.df[CH.df$temporal_order == "1", ]
  head(CH.df.tmp)
  list.patient_id.grt.three <- names(sort(table(CH.df.tmp$patient_id), decreasing = TRUE))[
    which( sort(table(CH.df.tmp$patient_id), decreasing = TRUE) >= 1)]
  rm(CH.df.tmp)

  CH.df.tmp <- CH.df
  vec.inc <- NULL
  for (i in 1: dim(CH.df.tmp)[1]){
    for (j in 1:length(list.patient_id.grt.three)){
      if (CH.df.tmp$patient_id[i] == list.patient_id.grt.three[j] ){
        vec.inc <- append(vec.inc, i)
      }
    }
  }
  CH.df.tmp <- CH.df.tmp[(vec.inc),]
  dim(CH.df.tmp)

  Nef.list <- Neff_Pr_List[[1]]$Nef.list
  Ter.list <- Neff_Pr_List[[1]]$Ter.list
  id.list.org <- Neff_Pr_List[[1]]$id.list.org
  id.list <- Neff_Pr_List[[1]]$id.list
  timePntday <- Neff_Pr_List[[1]]$Time.list
  timePntmon <- Neff_Pr_List[[1]]$TimeMonth.list

  mutCount.list <- Neff_Pr_List[[1]]$mutCount.list

  #pr.df.list <- vector (mode = "list", length = length(Pr.listS) )

  pr.df.list <- data.frame( Genes=c() , Func=c(), GeneFunc = c(), diff = c(), diffList=c(), diff2 = c(),
                            Pr=c(), Nef=c(), Treatment=c(), PtID=c(), day=c(),
                            month=c(), mutCount = c(), VAFdir=c(), PtID_org=c(), orderS=c()  )


  for (i in 1: length(Pr.listS)){ ## added for

    Nef.list[i]
    Ter.list[i]
    id.list[i]
    id.list[i]
    mutCount.list[i]
    CH.df.tmp2 <- CH.df.tmp[CH.df.tmp$patient_id == id.list.org[i],]
    CH.df.tmp2 <- CH.df.tmp2[order(CH.df.tmp2$amino_acid_change),]

    Pr.listS[[i]]

    diff <- rep(1, length(Pr.listS[[i]]))
    diff2 <- rep(1, length(Pr.listS[[i]]))

    diffList <- rep(1, length(Pr.listS[[i]]))

    VAFdir <- rep(1, length(Pr.listS[[i]]))

    pr.df <- data.frame(c( CH.df.tmp2$gene [CH.df.tmp2$temporal_order == "1"], "WT") ,
                        c(CH.df.tmp2$amino_acid_change [CH.df.tmp2$temporal_order == "1"], "WT"),
                        paste(c( CH.df.tmp2$gene [CH.df.tmp2$temporal_order == "1"], "WT") ,
                              c(CH.df.tmp2$amino_acid_change [CH.df.tmp2$temporal_order == "1"], "WT"),sep ="_" ),
                        Pr.listS[[i]], diff, diffList, diff2,
                        rep(Nef.list[i], length(Pr.listS[[i]])) ,
                        rep(Ter.list[i], length(Pr.listS[[i]])),
                        rep(id.list[i], length(Pr.listS[[i]])),
                        rep(timePntday[i], length(Pr.listS[[i]])),
                        rep(timePntmon[i], length(Pr.listS[[i]])),
                        rep(mutCount.list[i], length(Pr.listS[[i]])),
                        VAFdir,
                        rep( id.list.org[i], length(Pr.listS[[i]])),
                        rep( strsplit( id.list[i], "_")[[1]][3], length(Pr.listS[[i]]))
    )



    colnames(pr.df) <- c("Genes" , "Func", "GeneFunc",
                         "Pr",  "diff", "ground", "diff2",
                         "Nef", "Treatment", "PtID", "day",
                         "month", "mutCount",
                         "VAFdir", "PtID_org", "orderS" )

    orderStat <- strsplit( id.list[i], "_")[[1]][3]
    lsGene <- CH.df.tmp2$gene [CH.df.tmp2$temporal_order == 1 ]
    lsvar <- CH.df.tmp2$amino_acid_change [CH.df.tmp2$temporal_order == 1]

    for (p in 1:( dim(pr.df)[1] -1  ) ){

      CH.df.tmp2.gene <- CH.df.tmp2[CH.df.tmp2$gene ==lsGene[p], ]
      CH.df.tmp2.gene <- CH.df.tmp2.gene[CH.df.tmp2.gene$amino_acid_change ==lsvar[p], ]

      pr.df$VAFdir[p] <- as.numeric(CH.df.tmp2.gene$vaf [ which(CH.df.tmp2.gene$temporal_order == orderStat)  ]) -
        as.numeric(CH.df.tmp2.gene$vaf [ which(CH.df.tmp2.gene$temporal_order == 1) ] )
    }

    pr.df <- pr.df[order(pr.df$Pr, decreasing = F),]
    for (p in 1:( dim(pr.df)[1] -1  ) ){

      #pr.df$diff[p]  <- pr.df$Pr[p+1]/pr.df$Pr[p]
      pr.df$diff[p]  <- pr.df$Pr[dim(pr.df)[1]]/pr.df$Pr[p]
      pr.df$diff2[p]  <- pr.df$Pr[2] /pr.df$Pr[p]

      #pr.df$diffList[p] <- pr.df$GeneFunc[p+1]
      pr.df$ground[p] <- pr.df$GeneFunc[dim(pr.df)[1]]
    }

    #pr.df.list[i] <- pr.df
    pr.df.list <- bind_rows(pr.df.list, pr.df)

  } ## end for -- added

  prList.df <- pr.df.list
  prList.df$log10Diff <- log10(as.numeric(prList.df$diff))


  prList.df$dynamic <- NULL

  prList.df$dynamic [ which(prList.df$log10Diff >= SigOrderOfMagnitude_threshold) ] <- "sigPossible"
  prList.df$dynamic [ which(prList.df$log10Diff < SigOrderOfMagnitude_threshold) ] <- "noSelection"

  prList.df$dynamic [ intersect ( which (prList.df$dynamic == "sigPossible"), which (prList.df$VAFdir >= 0 ) ) ]  <- "Positive"
  prList.df$dynamic [ intersect ( which (prList.df$dynamic == "sigPossible"), which (prList.df$VAFdir <= 0 ) ) ]  <- "Negative"

  return (prList.df)
}
