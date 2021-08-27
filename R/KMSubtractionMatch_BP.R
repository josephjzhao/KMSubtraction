#' KMSubtractionMatch_BP
#'
#' This function conducts matching of patients between subgroup and an overall cohort through  minimum cost bipartite matching with the Hungarian algorithm.
#' @param df_overall Dataframe of the overall cohort containing time-to-event data. Prepared as "time" and "status" respectively.
#' @param df_subgroup Dataframe of the subgroup cohort containing time-to-event data. Prepared as "time" and "status" respectively.
#' @keywords 
#' @export
#' @examples
#' KMSubtractionMatch_BP(df_overall, df_subgroup)

KMSubtractionMatch_BP=function(df_overall, df_subgroup){
  #events
  overall=subset(df_overall, df_overall$status==1 )$time
  subgroup=subset(df_subgroup, df_subgroup$status==1 )$time
  matrix.events=sapply(subgroup, function(x){x-overall}, simplify=T) %>% abs()
  match.events=HungarianSolver(matrix.events)
  match.events$pairs[,2][match.events$pairs[,2]==0]=NA
  
  #censorships
  overall=subset(df_overall, df_overall$status==0 )$time
  subgroup=subset(df_subgroup, df_subgroup$status==0 )$time
  matrix.censorships=sapply(subgroup, function(x){x-overall}, simplify=T) %>% abs()
  match.censorships=HungarianSolver(matrix.censorships)
  match.censorships$pairs[,2][match.censorships$pairs[,2]==0]=NA
  
  df_match=rbind(
    rbind(cbind(subset(df_overall, df_overall$status==1), matchid=match.events$pairs[,2], curve="overall"), cbind(subset(df_subgroup, df_subgroup$status==1), matchid=1, curve="subgroup")),# events
    rbind(cbind(subset(df_overall, df_overall$status==0), matchid=match.censorships$pairs[,2], curve="overall"), cbind(subset(df_subgroup, df_subgroup$status==0), matchid=1, curve="subgroup")) # censors
  )
  
  df_match$subgroup=ifelse(is.na(df_match$matchid), 0, 1)
  df_match
}
