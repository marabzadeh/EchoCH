#' NeffAnalysis
#'
#' Calculate Neff for each patient with the Fitness likelihood for each variant in that patient
#'
#' @param CH.df a data.frame with CH mutations with the following columns::
#'                                              patient_id, gene, treatment, amino_acid_change,
#'                                              depth, vaf, temporal_order, pre/post_treatment,
#'                                             days_from_treatment_start_to_sample_collection )
#' @param mutNumber_threshold a threshhold on the number of mutations for the individual to be considered in the analysis (default 1)
#' @return Neff_Pr_List Neff_Pr_List[[1]]: a dataframe containing Neff results for each patient;
#'                      Neff_Pr_List[[2]]: a dataframe containing Pr results for each patient
#' @examples
#' Neff_Pr_List <- NeffAnalysis(CH.df)


NeffAnalysis <- function (CH.df, mutNumber_threshold = 1){

  CH.df$timepoint <- CH.df$days_from_treatment_start_to_sample_collection

  CH.df.tmp <- CH.df[CH.df$temporal_order == "1", ]
  list.ids.grt.three <- names(sort(table(CH.df.tmp$ids), decreasing = TRUE))[
   which( sort(table(CH.df.tmp$ids), decreasing = TRUE) >= mutNumber_threshold)]
  rm(CH.df.tmp)

  CH.df.tmp <- CH.df
  vec.inc <- NULL
  for (i in 1: dim(CH.df.tmp)[1]){
    for (j in 1:length(list.ids.grt.three)){
      if (CH.df.tmp$patient_id[i] == list.ids.grt.three[j] ){
        vec.inc <- append(vec.inc, i)
      }
    }
  }
  CH.df.tmp <- CH.df.tmp[(vec.inc),]
  dim(CH.df.tmp)

  ## For each go

  Nef.list <- NULL
  Time.list <- NULL
  Ter.list <- NULL
  Pr.listS <- NULL
  TimeMonth.list <- NULL
  order.list <- NULL
  mutCount.list <- NULL
  AF.list.pre <- NULL
  AF.list.post <- NULL
  GeneFunc <- NULL

  id.list <- NULL
  id.list.org <- NULL
  #var.list <- NULL

  count <- 1

  for (i in 1: length(list.ids.grt.three)) {

    CH.df.tmp2 <- CH.df.tmp[CH.df.tmp$patient_id == list.ids.grt.three[i],]
    #Ter.list[i] <- unique( CH.df.tmp2$treatment)

    s <- length(unique( CH.df.tmp2$amino_acid_change) ) + 1
    #mvec <- vector(mode = "integer", length = s)
    #nvec <- vector(mode = "integer", length = s)

    loops <- as.numeric(names(table(CH.df.tmp2$temporal_order)))

    if(length(loops) >= 2) {

      for (ll in 2: length(loops-1)) {

        startdate <- as.numeric(CH.df.tmp2$timepoint [ which(CH.df.tmp2$temporal_order == "1") ] ) [1]
        enddate <- as.numeric(CH.df.tmp2$timepoint [ which( CH.df.tmp2$temporal_order == loops[ll]) ] ) [1]

        CH.df.tmp2 <- CH.df.tmp2[order(CH.df.tmp2$amino_acid_change, CH.df.tmp2$temporal_order),]

        nvec <- as.numeric(CH.df.tmp2$vaf[which(CH.df.tmp2$temporal_order == "1") ] )
        mvec <- as.numeric(CH.df.tmp2$vaf[which(CH.df.tmp2$temporal_order == loops[ll]) ] )

        nvec <- append(nvec, 1- sum(nvec) )
        mvec <- append(mvec, 1- sum(mvec) )


        Nef <- calcN(s, mvec, nvec)
        #Nef

        Nef.list <- append(Nef.list, Nef)

        if (  Nef > 0){
          Pr.list <- vector(mode = "integer", length = s)
          for(j in 1:s){
            Pr.list[j] <- calcP(Nef,mvec[j],nvec[j])
          }
          cc <- length(Nef.list)
          Pr.listS[[cc]] <- Pr.list
        }



        Time.list <- append(Time.list, enddate - startdate)
        Ter.list <- append(Ter.list, unique( CH.df.tmp2$treatment))
        id.list <- append(id.list , paste0(list.ids.grt.three[i], "_", ll ) )
        id.list.org <- append(id.list.org , list.ids.grt.three[i] )
        order.list <- append(order.list , ll  )
        mutCount.list <- append(mutCount.list, (s-1) )

        AF.list.pre <- append(AF.list.pre, paste(round(nvec,5) ,collapse = ",") )
        AF.list.post <- append(AF.list.post, paste(round(mvec,5), collapse=",") )

        CH.df.tmp2$GeneFunc <- paste(c( CH.df.tmp2$GENE) ,
                                       c(CH.df.tmp2$amino_acid_change ),sep ="_" )

        GeneFunc.list <- c(CH.df.tmp2$GeneFunc[CH.df.tmp2$temporal_order == "1"], "WT")

        GeneFunc <- append(GeneFunc, paste(GeneFunc.list, collapse = ",") )


      }## end for loop
    }## end if loop
    count <- count +1
  }


  Nef.list <- unlist(Nef.list)
  Time.list <- unlist(Time.list)
  Ter.list <- unlist(Ter.list)
  id.list <- unlist(id.list)
  id.list.org <- unlist(id.list.org)
  order.list <- unlist(order.list)
  mutCount.list <- unlist(mutCount.list)
  AF.list.pre <- unlist(AF.list.pre)
  AF.list.post <- unlist(AF.list.post)
  GeneFunc <- unlist(GeneFunc)



  TimeMonth.list <- vector(mode = "character", length = length(Nef.list))

  for (i in 1: length(Nef.list)){

    if(Time.list[i] <= (18*30) ){
      TimeMonth.list[i] <- "<=18"
    } else if ( Time.list[i] > (18*30) ) {
      TimeMonth.list[i] <- ">18"
    }

  }


  df.Neff <- data.frame(Nef.list, Time.list, Ter.list, id.list, id.list.org, order.list,
                        TimeMonth.list,
                        mutCount.list,
                        AF.list.pre, AF.list.post, GeneFunc)

  ## ----------
  Neff_Pr_List <- vector (mode = "list", length = 2)
  Neff_Pr_List [[1]] <- df.Neff
  Neff_Pr_List [[2]] <- Pr.listS

  return (Neff_Pr_List)
}
