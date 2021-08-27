#' KMSubtractionMatch_BP
#'
#' This function conducts matching of patients between subgroup and an overall cohort through nearest neighbour matching with distances determined by logistic regression.
#' @param df_overall Dataframe of the overall cohort containing time-to-event data. Prepared as "time" and "status" respectively.
#' @param df_subgroup Dataframe of the subgroup cohort containing time-to-event data. Prepared as "time" and "status" respectively.
#' @keywords 
#' @export
#' @examples
#' KMSubtractionMatch_BP(df_overall, df_subgroup)

KMSubtractionMatch_Logit=function(df_overall, df_subgroup){
  
  df_overall$curve=0
  df_subgroup$curve=1
  
  df=bind_rows(df_overall, df_subgroup)
  
  match.events=matchit(curve~time, method="nearest", data=subset(df, df$status==1), ratio=1) 
  match.censorships=matchit(curve~time, method="nearest", data=subset(df, df$status==0), ratio=1)
  
  df_match=rbind(cbind(subset(df, df$status==1), matchid=match.events$subclass), 
                 cbind(subset(df, df$status==0), matchid=match.censorships$subclass))
  
  df_match$subgroup=ifelse(is.na(df_match$matchid), 0, 1)
  
  df_match$curve=ifelse(df_match$curve==0, "overall", "subgroup")
  
  df_match
  
}